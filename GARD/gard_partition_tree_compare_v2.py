#!/usr/bin/env python3

import os
import csv
import argparse
import subprocess
import copy

import matplotlib.pyplot as plt
from Bio import Phylo

# ----------------------------------------------------------------------------
# GARD partition tree comparison with gappy partition handling, sample renaming,
# plot pruning, IQ-TREE inference, midpoint-rooted tree plotting, and RF thresholds.
#
# Key behaviours:
#   - Splits whole-genome alignment at chosen GARD breakpoints.
#   - Removes all-missing / high-missingness taxa per partition before IQ-TREE.
#   - Renames selected taxa in FASTA files, IQ-TREE outputs, Newick files, and figures.
#   - Optionally removes selected taxa from plotted trees only (default: NC_005881.2).
#   - Runs IQ-TREE as `iqtree` by default, automatic ModelFinder model selection
#     (-m MFP), and 1000 ultrafast bootstraps by default.
#   - Midpoint-roots trees for output/plots; no outgroup rooting.
#   - Computes RF and normalized RF distances vs the per-partition pruned baseline.
#   - Adds threshold columns for normalized RF at 0.5, 0.8, and 0.95 by default.
#   - Generates IQ-TREE AU/SH topology-test commands.
# ----------------------------------------------------------------------------

MISSING_CHARS = set(['N', 'n', '-', '?'])

DEFAULT_RENAME_MAP = {
    '10-04648-0005': 'OsHV_NSW_2010_1',
    '12-03297-0004': 'OsHV_NSW_2011_1',
    '13-00205-0004': 'OsHV_NSW_2013_1',
    '19-00679-0186': 'OsHV_SA_2018_1',
    '16-00320-0005': 'OsHV_TAS_2016_1',
    '21-01515-0022': 'OsHV_TAS_2016_2',
    '21-01514-0003': 'OsHV_TAS_2018_1',
    '19-00509': 'OsHV_TAS_2019_1',
    '21-04170-0001': 'OsHV_TAS_2019_2',
    '24-00653-0004': 'OsHV_TAS_2024_1',
    '24-00653-0006': 'OsHV_TAS_2024_2',
    '24-00653-0013': 'OsHV_TAS_2024_3',
}


def parse_csv_mapping(path):
    """Read two-column TSV/CSV mapping: old_name <tab/comma> new_name."""
    mapping = {}
    if not path:
        return mapping
    with open(path, newline='') as fh:
        sample = fh.read(2048)
        fh.seek(0)
        dialect = csv.Sniffer().sniff(sample, delimiters='\t,')
        reader = csv.reader(fh, dialect)
        for row in reader:
            if not row or len(row) < 2:
                continue
            old = row[0].strip()
            new = row[1].strip()
            if not old or not new:
                continue
            if old.lower() in {'old', 'current', 'from', 'sample'}:
                continue
            mapping[old] = new
    return mapping


def rename_label(label, rename_map):
    return rename_map.get(label, label)


def rename_records(records, rename_map):
    out = []
    seen = set()
    for name, seq in records:
        new_name = rename_label(name, rename_map)
        if new_name in seen:
            raise RuntimeError(f'Duplicate sequence name after renaming: {new_name}')
        seen.add(new_name)
        out.append((new_name, seq))
    return out


def read_fasta(path):
    records = []
    header = None
    chunks = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip('\n')
            if not line:
                continue
            if line.startswith('>'):
                if header is not None:
                    records.append((header, ''.join(chunks)))
                header = line[1:].strip().split()[0]
                chunks = []
            else:
                chunks.append(line.strip())
    if header is not None:
        records.append((header, ''.join(chunks)))
    if not records:
        raise RuntimeError(f'No FASTA records found in {path}')
    lengths = {len(seq) for _, seq in records}
    if len(lengths) != 1:
        raise RuntimeError('Alignment sequences are not all the same length')
    return records, lengths.pop()


def write_fasta(records, path):
    with open(path, 'w') as out:
        for name, seq in records:
            out.write(f'>{name}\n')
            for i in range(0, len(seq), 80):
                out.write(seq[i:i+80] + '\n')


def parse_breakpoints(s):
    vals = []
    for x in s.split(','):
        x = x.strip()
        if x:
            vals.append(int(x))
    return sorted(set(vals))


def parse_thresholds(s):
    vals = []
    for x in s.split(','):
        x = x.strip()
        if x:
            vals.append(float(x))
    return vals


def split_alignment(records, aln_len, breakpoints):
    coords = []
    prev = 1
    for bp in breakpoints:
        coords.append((prev, bp))
        prev = bp + 1
    coords.append((prev, aln_len))
    partitions = []
    for i, (s, e) in enumerate(coords, start=1):
        part_records = [(name, seq[s-1:e]) for name, seq in records]
        partitions.append((i, s, e, part_records))
    return partitions


def missing_fraction(seq):
    if not seq:
        return 1.0
    return sum(1 for c in seq if c in MISSING_CHARS) / len(seq)


def informative_sites(records):
    if not records:
        return 0, 0
    L = len(records[0][1])
    if any(len(seq) != L for _, seq in records):
        raise RuntimeError('Partition sequences not same length')
    usable = 0
    variable = 0
    for i in range(L):
        resolved = [seq[i].upper() for _, seq in records if seq[i] not in MISSING_CHARS]
        if len(resolved) >= 2:
            usable += 1
            if len(set(resolved)) >= 2:
                variable += 1
    return usable, variable


def filter_partition_records(records, drop_missing_above=1.0, remove_from_analysis=None):
    remove_from_analysis = set(remove_from_analysis or [])
    kept = []
    removed = []
    for name, seq in records:
        frac = missing_fraction(seq)
        if name in remove_from_analysis:
            removed.append((name, frac, 'explicitly_removed_from_analysis'))
        elif frac >= drop_missing_above:
            removed.append((name, frac, 'all_missing_or_above_threshold'))
        else:
            kept.append((name, seq))
    return kept, removed


def run_iqtree(aln_fasta, prefix, executable='iqtree', threads=1, model='MFP', ufboot=1000, extra_args=None):
    cmd = [executable, '-s', aln_fasta, '-pre', prefix, '-nt', str(threads), '-m', model, '-st', 'DNA']
    if ufboot and int(ufboot) > 0:
        cmd.extend(['-B', str(int(ufboot))])
    if extra_args:
        cmd.extend(extra_args)
    subprocess.run(cmd, check=True)
    treefile = prefix + '.treefile'
    if not os.path.isfile(treefile):
        raise RuntimeError(f'IQ-TREE finished but treefile not found: {treefile}')
    return treefile


def load_tree(path):
    return Phylo.read(path, 'newick')


def terminal_names(tree):
    return sorted(t.name for t in tree.get_terminals())


def rename_tree_tips(tree, rename_map):
    tree = copy.deepcopy(tree)
    seen = set()
    for term in tree.get_terminals():
        term.name = rename_label(term.name, rename_map)
        if term.name in seen:
            raise RuntimeError(f'Duplicate tree tip after renaming: {term.name}')
        seen.add(term.name)
    return tree


def midpoint_rooted_tree(in_tree):
    tree = copy.deepcopy(in_tree)
    try:
        tree.root_at_midpoint()
    except Exception:
        pass
    return tree


def write_tree(tree, out_path):
    Phylo.write(tree, out_path, 'newick')
    return out_path


def midpoint_root_file(tree_path, out_path, rename_map=None):
    tree = load_tree(tree_path)
    if rename_map:
        tree = rename_tree_tips(tree, rename_map)
    rooted = midpoint_rooted_tree(tree)
    return write_tree(rooted, out_path)


def prune_tree_to_taxa(tree, keep_taxa):
    tree = copy.deepcopy(tree)
    keep_taxa = set(keep_taxa)
    for term in list(tree.get_terminals()):
        if term.name not in keep_taxa:
            tree.prune(term)
    return tree


def prune_tree_for_plot(tree, remove_taxa):
    tree = copy.deepcopy(tree)
    remove_taxa = set(remove_taxa or [])
    for term in list(tree.get_terminals()):
        if term.name in remove_taxa:
            try:
                tree.prune(term)
            except Exception:
                pass
    return tree


def bipartitions(tree, shared_leaves=None):
    leaves_all = set(terminal_names(tree))
    if shared_leaves is None:
        shared_leaves = leaves_all
    shared_leaves = set(shared_leaves)
    root = tree.root
    parts = set()
    total = frozenset(shared_leaves)

    def descendants(clade):
        return set(t.name for t in clade.get_terminals() if t.name in shared_leaves)

    for clade in tree.find_clades(order='preorder'):
        if clade == root:
            continue
        desc = frozenset(descendants(clade))
        if len(desc) == 0 or len(desc) == len(total):
            continue
        other = total - desc
        canonical = desc if len(desc) < len(other) else other
        if 1 < len(canonical) < len(total) - 1:
            parts.add(canonical)
    return parts



def parse_node_support(clade, paired_mode='min'):
    """Read IQ-TREE support as a 0..1 value from UFBoot or SH-aLRT/UFBoot."""
    import re
    raw = clade.name if clade.name not in (None, '') else clade.confidence
    if raw is None:
        return None
    vals = [float(x) for x in re.findall(r'\d+(?:\.\d+)?', str(raw))]
    if not vals:
        return None
    vals = [v / 100.0 if v > 1 else v for v in vals]
    if len(vals) == 1:
        return vals[0]
    if paired_mode == 'ufboot':
        return vals[-1]
    if paired_mode == 'sh-alrt':
        return vals[0]
    if paired_mode == 'mean':
        return sum(vals[:2]) / 2.0
    return min(vals[:2])


def supported_bipartitions(tree, threshold, paired_mode='min'):
    """Return non-trivial unrooted splits whose node support meets threshold."""
    leaves = set(terminal_names(tree))
    total = frozenset(leaves)
    parts = set()
    for clade in tree.find_clades(order='preorder'):
        if clade == tree.root or clade.is_terminal():
            continue
        support = parse_node_support(clade, paired_mode)
        if support is None or support < threshold:
            continue
        desc = frozenset(t.name for t in clade.get_terminals())
        other = total - desc
        if len(desc) < 2 or len(other) < 2:
            continue
        if len(desc) < len(other):
            canonical = desc
        elif len(other) < len(desc):
            canonical = other
        else:
            canonical = min((desc, other), key=lambda x: tuple(sorted(x)))
        parts.add(canonical)
    return parts


def compare_split_sets(a, b):
    shared = len(a & b)
    union = len(a | b)
    rf = len(a - b) + len(b - a)
    denom = len(a) + len(b)
    return {
        'rf': rf,
        'normalized_rf': rf / denom if denom else 0.0,
        'shared_splits': shared,
        'baseline_splits': len(a),
        'partition_splits': len(b),
        'jaccard': shared / union if union else 1.0,
        'precision': shared / len(b) if b else (1.0 if not a else 0.0),
        'recall': shared / len(a) if a else (1.0 if not b else 0.0),
    }


def support_threshold_comparisons(tree1, tree2, thresholds, paired_mode='min'):
    """Prune BOTH trees to identical taxa, then compare supported splits."""
    shared = set(terminal_names(tree1)) & set(terminal_names(tree2))
    if len(shared) < 4:
        return None, None, len(shared), {}
    t1 = prune_tree_to_taxa(tree1, shared)
    t2 = prune_tree_to_taxa(tree2, shared)
    if set(terminal_names(t1)) != set(terminal_names(t2)):
        raise RuntimeError('Taxon sets differ after matched pruning')
    output = {}
    for threshold in thresholds:
        output[threshold] = compare_split_sets(
            supported_bipartitions(t1, threshold, paired_mode),
            supported_bipartitions(t2, threshold, paired_mode),
        )
    return t1, t2, len(shared), output

def rf_distance(tree1, tree2):
    leaves1 = set(terminal_names(tree1))
    leaves2 = set(terminal_names(tree2))
    shared = leaves1 & leaves2
    if len(shared) < 4:
        return None, None, len(shared)
    b1 = bipartitions(tree1, shared_leaves=shared)
    b2 = bipartitions(tree2, shared_leaves=shared)
    rf = len(b1 - b2) + len(b2 - b1)
    max_rf = len(b1) + len(b2)
    norm = rf / max_rf if max_rf else 0.0
    return rf, norm, len(shared)


def rf_threshold_class(norm_rf, thresholds):
    if norm_rf in ('', None):
        return ''
    passed = [t for t in thresholds if norm_rf >= t]
    if not passed:
        return f'<{thresholds[0]}' if thresholds else ''
    return f'>={max(passed)}'


def draw_tree_to_png(tree_path, out_png, title=None, midpoint_root=True, remove_from_plot=None, rename_map=None):
    tree = load_tree(tree_path)
    if rename_map:
        tree = rename_tree_tips(tree, rename_map)
    if remove_from_plot:
        tree = prune_tree_for_plot(tree, remove_from_plot)
    if midpoint_root:
        tree = midpoint_rooted_tree(tree)
    tree.ladderize()
    n_tips = len(tree.get_terminals())
    height = max(8, 0.28 * n_tips)
    fig = plt.figure(figsize=(12, height))
    ax = fig.add_subplot(1, 1, 1)
    Phylo.draw(tree, do_show=False, axes=ax)
    if title:
        ax.set_title(title)
    fig.savefig(out_png, dpi=220, bbox_inches='tight')
    plt.close(fig)


def write_iqtree_test_commands(partitions, out_sh, executable='iqtree'):
    with open(out_sh, 'w') as out:
        out.write('#!/usr/bin/env bash\nset -euo pipefail\n\n')
        out.write('# IQ-TREE topology tests (AU/SH family) comparing each partition tree to the partition-pruned baseline tree.\n')
        out.write('MODEL=${MODEL:-MFP}\n')
        for part in partitions:
            aln = part.get('alignment')
            tree = part.get('treefile')
            baseline_pruned = part.get('baseline_pruned_tree')
            if not aln or not tree or not baseline_pruned:
                continue
            prefix = part['prefix'] + '.topotest'
            treeset = part['prefix'] + '.candidate_trees.nwk'
            out.write(f'cat {baseline_pruned} {tree} > {treeset}\n')
            out.write(f'{executable} -s {aln} -st DNA -m $MODEL -z {treeset} -n 0 -zb 10000 -au -pre {prefix}\n\n')


def comma_list(s):
    if not s:
        return []
    return [x.strip() for x in s.split(',') if x.strip()]


def main():
    ap = argparse.ArgumentParser(description='Split an alignment at chosen GARD breakpoints, infer per-partition IQ-TREE trees, rename tips, hide selected taxa from plots, and compute RF threshold summaries.')
    ap.add_argument('--alignment', required=True, help='Whole-alignment FASTA used for GARD')
    ap.add_argument('--breakpoints', required=True, help='Comma-separated breakpoint list, e.g. 8940,42522,66255')
    ap.add_argument('--baseline-tree', required=True, help='Zero-breakpoint baseline tree (Newick) for the same sample set')
    ap.add_argument('--outdir', required=True, help='Output directory')
    ap.add_argument('--run-iqtree', action='store_true', help='Infer trees for each partition using IQ-TREE')
    ap.add_argument('--threads', type=int, default=1)
    ap.add_argument('--iqtree-exe', default='iqtree', help='IQ-TREE executable name/path (default: iqtree)')
    ap.add_argument('--iqtree-model', default='MFP', help='IQ-TREE model setting. Default MFP = automatic ModelFinder model selection.')
    ap.add_argument('--ufboot', type=int, default=1000, help='Ultrafast bootstrap replicates per tree (default 1000; set 0 to disable)')
    ap.add_argument('--iqtree-extra', default='', help='Extra raw arguments for IQ-TREE, e.g. "-bnni"')
    ap.add_argument('--drop-missing-above', type=float, default=1.0, help='Drop taxa from a partition if fraction of N,-,? is >= this value')
    ap.add_argument('--min-taxa', type=int, default=4, help='Skip partition if fewer than this many taxa remain after filtering')
    ap.add_argument('--min-usable-sites', type=int, default=20, help='Skip partition if fewer than this many usable columns remain')
    ap.add_argument('--rf-thresholds', default='0.5,0.8,0.95', help='Comma-separated normalized RF category thresholds')
    ap.add_argument('--support-thresholds', default='0.5,0.8,0.9', help='Node-support thresholds for supported-split RF comparisons')
    ap.add_argument('--paired-support-mode', choices=['min','ufboot','sh-alrt','mean'], default='min', help='For SH-aLRT/UFBoot labels, min requires both tests to meet threshold')
    ap.add_argument('--rename-map', default='', help='Optional TSV/CSV with old_name,new_name. Built-in OsHV rename map is used by default.')
    ap.add_argument('--no-default-renames', action='store_true', help='Disable built-in OsHV sample renaming')
    ap.add_argument('--remove-from-plots', default='NC_005881.2', help='Comma-separated taxa to remove from plotted PNG trees only')
    ap.add_argument('--remove-from-analysis', default='', help='Comma-separated taxa to remove from partition FASTAs / IQ-TREE / RF comparisons')
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    for subdir in ['alignments', 'trees', 'plots', 'logs']:
        os.makedirs(os.path.join(args.outdir, subdir), exist_ok=True)

    rename_map = {} if args.no_default_renames else dict(DEFAULT_RENAME_MAP)
    rename_map.update(parse_csv_mapping(args.rename_map))
    remove_from_plots = {rename_label(x, rename_map) for x in comma_list(args.remove_from_plots)}
    remove_from_analysis = {rename_label(x, rename_map) for x in comma_list(args.remove_from_analysis)}
    thresholds = parse_thresholds(args.rf_thresholds)
    support_thresholds = parse_thresholds(args.support_thresholds)
    extra_args = args.iqtree_extra.split() if args.iqtree_extra.strip() else None

    original_records, aln_len = read_fasta(args.alignment)
    records = rename_records(original_records, rename_map)
    breakpoints = parse_breakpoints(args.breakpoints)
    partitions = split_alignment(records, aln_len, breakpoints)

    # Rename baseline once and use this internally for all comparisons.
    baseline_tree_orig = load_tree(args.baseline_tree)
    baseline_tree_renamed = rename_tree_tips(baseline_tree_orig, rename_map)
    baseline_renamed_path = os.path.join(args.outdir, 'trees', 'baseline_tree.renamed.nwk')
    write_tree(baseline_tree_renamed, baseline_renamed_path)

    baseline_mid = os.path.join(args.outdir, 'trees', 'baseline_tree.renamed.midpoint_rooted.nwk')
    write_tree(midpoint_rooted_tree(baseline_tree_renamed), baseline_mid)
    draw_tree_to_png(baseline_renamed_path, os.path.join(args.outdir, 'plots', 'baseline_tree.png'), title='Baseline (0-breakpoint) tree', midpoint_root=True, remove_from_plot=remove_from_plots)

    summary_rows = []
    removed_rows = []
    part_meta = []

    for idx, start, end, part_records in partitions:
        part_name = f'partition_{idx}_{start}_{end}'
        filtered_records, removed = filter_partition_records(part_records, drop_missing_above=args.drop_missing_above, remove_from_analysis=remove_from_analysis)
        for name, frac, reason in removed:
            removed_rows.append({'partition': idx, 'start': start, 'end': end, 'taxon': name, 'missing_fraction': frac, 'reason': reason})

        usable_sites, variable_sites = informative_sites(filtered_records) if filtered_records else (0, 0)
        keep_taxa = [name for name, _ in filtered_records]

        baseline_pruned = os.path.join(args.outdir, 'trees', part_name + '.baseline_pruned.renamed.nwk')
        baseline_pruned_mid = os.path.join(args.outdir, 'trees', part_name + '.baseline_pruned.renamed.midpoint_rooted.nwk')
        if keep_taxa:
            b_pruned = prune_tree_to_taxa(baseline_tree_renamed, keep_taxa)
            write_tree(b_pruned, baseline_pruned)
            write_tree(midpoint_rooted_tree(b_pruned), baseline_pruned_mid)
        else:
            baseline_pruned = ''
            baseline_pruned_mid = ''

        aln_path = os.path.join(args.outdir, 'alignments', part_name + '.renamed.fasta')
        if filtered_records:
            write_fasta(filtered_records, aln_path)
        else:
            aln_path = ''

        status = 'ready'
        details = ''
        if len(filtered_records) < args.min_taxa:
            status = 'skipped_too_few_taxa'
            details = f'{len(filtered_records)} taxa remain after filtering'
        elif usable_sites < args.min_usable_sites:
            status = 'skipped_too_few_usable_sites'
            details = f'{usable_sites} usable sites remain after filtering'

        treefile = ''
        treefile_mid = ''
        rf = ''
        norm_rf = ''
        n_shared = ''
        support_results = {}
        matched_base_path = ''
        matched_part_path = ''

        if status == 'ready' and args.run_iqtree:
            try:
                prefix = os.path.join(args.outdir, 'trees', part_name)
                treefile = run_iqtree(aln_path, prefix, executable=args.iqtree_exe, threads=args.threads, model=args.iqtree_model, ufboot=args.ufboot, extra_args=extra_args)
                # IQ-TREE read the renamed FASTA, so treefile labels should already be renamed.
                p_tree = load_tree(treefile)
                treefile_mid = os.path.join(args.outdir, 'trees', part_name + '.midpoint_rooted.nwk')
                write_tree(midpoint_rooted_tree(p_tree), treefile_mid)
                draw_tree_to_png(treefile, os.path.join(args.outdir, 'plots', part_name + '.png'), title=part_name, midpoint_root=True, remove_from_plot=remove_from_plots)

                b_tree = load_tree(baseline_pruned)
                # Match BOTH trees to the exact intersecting taxon set before RF.
                matched_base, matched_part, shared_val, support_results = support_threshold_comparisons(
                    b_tree, p_tree, support_thresholds, args.paired_support_mode)
                if matched_base is not None:
                    matched_base_path = os.path.join(args.outdir, 'trees', part_name + '.baseline_taxon_matched.nwk')
                    matched_part_path = os.path.join(args.outdir, 'trees', part_name + '.partition_taxon_matched.nwk')
                    write_tree(matched_base, matched_base_path)
                    write_tree(matched_part, matched_part_path)
                    rf_val, norm_val, shared_val = rf_distance(matched_base, matched_part)
                else:
                    matched_base_path = ''
                    matched_part_path = ''
                    rf_val, norm_val = None, None
                rf = '' if rf_val is None else rf_val
                norm_rf = '' if norm_val is None else norm_val
                n_shared = '' if shared_val is None else shared_val
            except subprocess.CalledProcessError as e:
                status = 'iqtree_failed'
                details = f'IQ-TREE failed with exit status {e.returncode}'

        row = {
            'partition': idx,
            'start': start,
            'end': end,
            'length': end - start + 1,
            'n_taxa_before_filtering': len(part_records),
            'n_taxa_after_filtering': len(filtered_records),
            'n_removed_for_missingness_or_analysis': len(removed),
            'usable_sites_after_filtering': usable_sites,
            'variable_sites_after_filtering': variable_sites,
            'status': status,
            'details': details,
            'treefile': treefile,
            'treefile_midpoint_rooted': treefile_mid,
            'baseline_pruned_tree': baseline_pruned,
            'baseline_pruned_midpoint_rooted': baseline_pruned_mid,
            'rf_distance_vs_baseline': rf,
            'normalized_rf_vs_baseline': norm_rf,
            'rf_threshold_class': rf_threshold_class(norm_rf, thresholds),
            'shared_taxa_with_baseline': n_shared,
            'baseline_taxon_matched_tree': matched_base_path,
            'partition_taxon_matched_tree': matched_part_path,
        }
        for threshold, metrics in support_results.items():
            safe_support = str(threshold).replace('.', '_')
            for metric, value in metrics.items():
                row[f'support_{safe_support}_{metric}'] = value
        for t in thresholds:
            safe = str(t).replace('.', '_')
            row[f'normalized_rf_ge_{safe}'] = 'yes' if norm_rf != '' and norm_rf >= t else ('no' if norm_rf != '' else '')
        summary_rows.append(row)

        part_meta.append({
            'partition': idx,
            'alignment': aln_path,
            'treefile': treefile,
            'baseline_pruned_tree': baseline_pruned,
            'prefix': os.path.join(args.outdir, 'trees', part_name),
        })

    summary_tsv = os.path.join(args.outdir, 'partition_tree_summary.tsv')
    threshold_fields = [f'normalized_rf_ge_{str(t).replace(".", "_")}' for t in thresholds]
    support_fields = []
    for t in support_thresholds:
        safe = str(t).replace('.', '_')
        for metric in ['rf','normalized_rf','shared_splits','baseline_splits','partition_splits','jaccard','precision','recall']:
            support_fields.append(f'support_{safe}_{metric}')
    fields = [
        'partition', 'start', 'end', 'length',
        'n_taxa_before_filtering', 'n_taxa_after_filtering', 'n_removed_for_missingness_or_analysis',
        'usable_sites_after_filtering', 'variable_sites_after_filtering',
        'status', 'details',
        'treefile', 'treefile_midpoint_rooted', 'baseline_pruned_tree', 'baseline_pruned_midpoint_rooted',
        'baseline_taxon_matched_tree', 'partition_taxon_matched_tree',
        'rf_distance_vs_baseline', 'normalized_rf_vs_baseline', 'rf_threshold_class'
    ] + threshold_fields + ['shared_taxa_with_baseline'] + support_fields

    with open(summary_tsv, 'w', newline='') as out:
        writer = csv.DictWriter(out, fieldnames=fields, delimiter='\t')
        writer.writeheader()
        for row in summary_rows:
            writer.writerow(row)

    removed_tsv = os.path.join(args.outdir, 'logs', 'removed_sequences_by_partition.tsv')
    with open(removed_tsv, 'w', newline='') as out:
        writer = csv.DictWriter(out, fieldnames=['partition', 'start', 'end', 'taxon', 'missing_fraction', 'reason'], delimiter='\t')
        writer.writeheader()
        for row in removed_rows:
            writer.writerow(row)

    map_tsv = os.path.join(args.outdir, 'logs', 'applied_tip_rename_map.tsv')
    with open(map_tsv, 'w', newline='') as out:
        writer = csv.writer(out, delimiter='\t')
        writer.writerow(['old_name', 'new_name'])
        for old, new in sorted(rename_map.items()):
            writer.writerow([old, new])

    plot_removed_tsv = os.path.join(args.outdir, 'logs', 'taxa_removed_from_plots_only.tsv')
    with open(plot_removed_tsv, 'w', newline='') as out:
        writer = csv.writer(out, delimiter='\t')
        writer.writerow(['taxon'])
        for taxon in sorted(remove_from_plots):
            writer.writerow([taxon])

    iqtree_cmds = os.path.join(args.outdir, 'run_iqtree_topology_tests.sh')
    write_iqtree_test_commands(part_meta, iqtree_cmds, executable=args.iqtree_exe)
    os.chmod(iqtree_cmds, 0o755)

    print('Done.')
    print(f'Partition summary: {summary_tsv}')
    print(f'Removed-sequence log: {removed_tsv}')
    print(f'Applied rename map: {map_tsv}')
    print(f'Taxa removed from plots only: {plot_removed_tsv}')
    print(f'Baseline renamed midpoint-rooted tree: {baseline_mid}')
    print(f'Baseline tree plot: {os.path.join(args.outdir, "plots", "baseline_tree.png")}')
    print(f'IQ-TREE topology-test command script: {iqtree_cmds}')


if __name__ == '__main__':
    main()
