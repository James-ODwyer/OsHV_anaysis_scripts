#!/usr/bin/env python3

import os
import re
import glob
import csv
import argparse
from collections import defaultdict

ORF_REGEX = re.compile(r'\bORF[_\-\s]?[A-Za-z0-9.]+\b', re.IGNORECASE)


def normalize_gene_name(name):
    if name is None:
        return None
    name = str(name).strip()
    if not name or name.upper() == 'NA':
        return None
    m = ORF_REGEX.search(name)
    if m:
        x = m.group(0).upper()
        x = x.replace('ORF_', 'ORF').replace('ORF-', 'ORF').replace(' ', '')
        return x
    return re.sub(r'\s+', '_', name)


def parse_attrs(attr_str):
    attrs = {}
    for item in attr_str.strip().split(';'):
        if '=' in item:
            k, v = item.split('=', 1)
            attrs[k] = v
    return attrs


def extract_gene_name(attrs):
    for key in ['gene', 'Name', 'locus_tag', 'product', 'note', 'ID']:
        if key in attrs:
            g = normalize_gene_name(attrs[key])
            if g:
                return g
    return None


def parse_run_arg(run_arg):
    if '=' not in run_arg:
        raise ValueError(f'Run argument must be LABEL=PATH, got: {run_arg}')
    label, path = run_arg.split('=', 1)
    label = label.strip(); path = path.strip()
    if not label or not path:
        raise ValueError(f'Invalid run argument: {run_arg}')
    return label, path


def resolve_primary_reference_gff(run_dir, label):
    ref_dir = os.path.join(run_dir, '01_reference', 'ref_gffs')
    if not os.path.isdir(ref_dir):
        raise FileNotFoundError(f'Missing reference GFF directory: {ref_dir}')
    exact = os.path.join(ref_dir, f'{label}.gff3')
    if os.path.isfile(exact):
        return exact
    matches = sorted(glob.glob(os.path.join(ref_dir, f'{label}*.gff3')))
    if len(matches) == 1:
        return matches[0]
    if len(matches) > 1:
        raise RuntimeError(f"Multiple reference GFFs match label '{label}': {matches}")
    raise FileNotFoundError(f"Could not find primary reference GFF for label '{label}' in {ref_dir}")


def compute_gene_lengths_from_gff(gff_file):
    gene_id_to_name = {}
    mrna_to_gene = {}
    cds_by_gene = defaultdict(int)
    with open(gff_file) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            toks = line.rstrip('\n').split('\t')
            if len(toks) != 9:
                continue
            seqid, source, ftype, start, end, score, strand, phase, attrs_str = toks
            start = int(start); end = int(end)
            attrs = parse_attrs(attrs_str)
            if ftype == 'gene':
                gid = attrs.get('ID'); gname = extract_gene_name(attrs)
                if gid and gname:
                    gene_id_to_name[gid] = gname
            elif ftype in {'mRNA', 'transcript'}:
                mid = attrs.get('ID'); parent = attrs.get('Parent')
                if mid and parent:
                    mrna_to_gene[mid] = parent
            elif ftype == 'CDS':
                parent = attrs.get('Parent')
                gname = None
                if parent in gene_id_to_name:
                    gname = gene_id_to_name[parent]
                elif parent in mrna_to_gene and mrna_to_gene[parent] in gene_id_to_name:
                    gname = gene_id_to_name[mrna_to_gene[parent]]
                else:
                    gname = extract_gene_name(attrs)
                if gname:
                    cds_by_gene[gname] += (end - start + 1)
    return {g: {'cds_nt': nt, 'protein_aa': nt // 3} for g, nt in cds_by_gene.items()}


def find_final_gffs(run_dir):
    pattern = os.path.join(run_dir, '05_final', '*', '*.final.gff3')
    return {os.path.basename(g).replace('.final.gff3', ''): g for g in sorted(glob.glob(pattern))}


def find_truncate_summary_files(run_dir):
    pattern = os.path.join(run_dir, '06_qc', '*', '*.truncate_aware.tsv')
    return {os.path.basename(f).replace('.truncate_aware.tsv', ''): f for f in sorted(glob.glob(pattern))}


def parse_truncate_summary(tsv_file):
    rows = []
    with open(tsv_file, newline='') as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        for row in reader:
            rows.append(row)
    return rows


def fmt_value_or_range(values):
    vals = sorted(set(v for v in values if v is not None))
    if not vals:
        return ''
    if len(vals) == 1:
        return str(vals[0])
    return f'{vals[0]}-{vals[-1]}'


def ensure_dir(path):
    os.makedirs(path, exist_ok=True)


def main():
    parser = argparse.ArgumentParser(description='Summarise genes with definitive premature stop codons, excluding cases where any reference-based run shows the same sample/gene as complete.')
    parser.add_argument('--run', action='append', required=True, help='Run specification LABEL=PATH; use multiple times')
    parser.add_argument('--outdir', required=True)
    parser.add_argument('--complete-fraction', type=float, default=0.95, help='Call a sample/gene complete if any run has observed CDS length >= this fraction of that run\'s reference CDS length (default: 0.95)')
    args = parser.parse_args()

    ensure_dir(args.outdir)
    runs = {}
    for ra in args.run:
        label, path = parse_run_arg(ra)
        if label in runs:
            raise ValueError(f'Duplicate label: {label}')
        if not os.path.isdir(path):
            raise FileNotFoundError(f'Run directory not found: {path}')
        runs[label] = path

    ref_lengths_by_run = {}
    final_lengths_by_run = {}
    for label, run_dir in runs.items():
        ref_gff = resolve_primary_reference_gff(run_dir, label)
        ref_lengths_by_run[label] = compute_gene_lengths_from_gff(ref_gff)
        final_lengths_by_run[label] = {sample: compute_gene_lengths_from_gff(gff) for sample, gff in find_final_gffs(run_dir).items()}

    # Determine which sample/gene pairs appear complete in ANY run.
    complete_in_any = set()
    for label, per_sample in final_lengths_by_run.items():
        for sample, gene_lengths in per_sample.items():
            for gene, obs in gene_lengths.items():
                ref_nt = ref_lengths_by_run[label].get(gene, {}).get('cds_nt')
                if ref_nt and obs['cds_nt'] >= ref_nt * args.complete_fraction:
                    complete_in_any.add((sample, gene))

    event_rows = []
    skipped_complete_rows = []
    agg = defaultdict(lambda: {
        'run_labels': set(),
        'ref_nts': [],
        'ref_aas': [],
        'obs_nts': [],
        'obs_aas': [],
        'trim_bps': [],
    })

    for label, run_dir in runs.items():
        final_gffs = find_final_gffs(run_dir)
        truncate_summaries = find_truncate_summary_files(run_dir)
        for sample, summary_file in sorted(truncate_summaries.items()):
            if sample not in final_gffs:
                continue
            final_lengths = final_lengths_by_run[label].get(sample, {})
            rows = parse_truncate_summary(summary_file)
            for row in rows:
                gene = normalize_gene_name(row.get('gene'))
                if not gene:
                    continue
                try:
                    trim_bp = int(row.get('internal_stop_trim_bp', 0) or 0)
                except ValueError:
                    trim_bp = 0
                if trim_bp <= 0:
                    continue

                # if any run shows same sample/gene complete, suppress this false positive
                if (sample, gene) in complete_in_any:
                    skipped_complete_rows.append({'run_label': label, 'sample': sample, 'gene': gene, 'reason': 'complete_in_at_least_one_reference_run'})
                    continue

                ref_nt = ref_lengths_by_run[label].get(gene, {}).get('cds_nt')
                ref_aa = ref_lengths_by_run[label].get(gene, {}).get('protein_aa')
                obs_nt = final_lengths.get(gene, {}).get('cds_nt')
                obs_aa = final_lengths.get(gene, {}).get('protein_aa')
                event_rows.append({
                    'run_label': label,
                    'sample': sample,
                    'gene': gene,
                    'reference_cds_nt': ref_nt,
                    'reference_protein_aa': ref_aa,
                    'observed_cds_nt': obs_nt,
                    'observed_protein_aa': obs_aa,
                    'internal_stop_trim_bp': trim_bp,
                })
                key = (sample, gene)
                agg[key]['run_labels'].add(label)
                if ref_nt is not None:
                    agg[key]['ref_nts'].append(ref_nt)
                if ref_aa is not None:
                    agg[key]['ref_aas'].append(ref_aa)
                if obs_nt is not None:
                    agg[key]['obs_nts'].append(obs_nt)
                if obs_aa is not None:
                    agg[key]['obs_aas'].append(obs_aa)
                agg[key]['trim_bps'].append(trim_bp)

    event_out = os.path.join(args.outdir, 'premature_stop_per_run_sample.tsv')
    with open(event_out, 'w', newline='') as out:
        writer = csv.DictWriter(out, fieldnames=['run_label', 'sample', 'gene', 'reference_cds_nt', 'reference_protein_aa', 'observed_cds_nt', 'observed_protein_aa', 'internal_stop_trim_bp'], delimiter='\t')
        writer.writeheader()
        for row in sorted(event_rows, key=lambda x: (x['gene'], x['sample'], x['run_label'])):
            writer.writerow(row)

    collapsed_out = os.path.join(args.outdir, 'premature_stop_unique_sample_gene.tsv')
    gene_summary_out = os.path.join(args.outdir, 'premature_stop_gene_summary.tsv')
    suppressed_out = os.path.join(args.outdir, 'suppressed_complete_in_any_reference.tsv')

    gene_summary = defaultdict(lambda: {'samples': set(), 'run_labels': set(), 'ref_nts': [], 'ref_aas': [], 'obs_nts': [], 'obs_aas': [], 'sample_details': []})

    with open(collapsed_out, 'w', newline='') as out:
        writer = csv.DictWriter(out, fieldnames=['gene', 'sample', 'runs_detected', 'reference_cds_nt', 'reference_protein_aa', 'observed_cds_nt', 'observed_protein_aa', 'internal_stop_trim_bp'], delimiter='\t')
        writer.writeheader()
        for (sample, gene), d in sorted(agg.items(), key=lambda x: (x[0][1], x[0][0])):
            row = {
                'gene': gene,
                'sample': sample,
                'runs_detected': ','.join(sorted(d['run_labels'])),
                'reference_cds_nt': fmt_value_or_range(d['ref_nts']),
                'reference_protein_aa': fmt_value_or_range(d['ref_aas']),
                'observed_cds_nt': fmt_value_or_range(d['obs_nts']),
                'observed_protein_aa': fmt_value_or_range(d['obs_aas']),
                'internal_stop_trim_bp': fmt_value_or_range(d['trim_bps'])
            }
            writer.writerow(row)
            gs = gene_summary[gene]
            gs['samples'].add(sample)
            gs['run_labels'].update(d['run_labels'])
            gs['ref_nts'].extend(d['ref_nts'])
            gs['ref_aas'].extend(d['ref_aas'])
            gs['obs_nts'].extend(d['obs_nts'])
            gs['obs_aas'].extend(d['obs_aas'])
            gs['sample_details'].append(f"{sample}:{fmt_value_or_range(d['obs_aas'])}aa/{fmt_value_or_range(d['obs_nts'])}nt")

    with open(gene_summary_out, 'w', newline='') as out:
        writer = csv.writer(out, delimiter='\t')
        writer.writerow(['gene', 'n_samples_truncated', 'samples', 'runs_detected', 'reference_cds_nt', 'reference_protein_aa', 'observed_cds_nt', 'observed_protein_aa', 'sample_observed_lengths'])
        for gene, d in sorted(gene_summary.items()):
            writer.writerow([
                gene,
                len(d['samples']),
                ','.join(sorted(d['samples'])),
                ','.join(sorted(d['run_labels'])),
                fmt_value_or_range(d['ref_nts']),
                fmt_value_or_range(d['ref_aas']),
                fmt_value_or_range(d['obs_nts']),
                fmt_value_or_range(d['obs_aas']),
                ';'.join(sorted(d['sample_details']))
            ])

    with open(suppressed_out, 'w', newline='') as out:
        writer = csv.DictWriter(out, fieldnames=['run_label', 'sample', 'gene', 'reason'], delimiter='\t')
        writer.writeheader()
        for row in sorted(skipped_complete_rows, key=lambda x: (x['gene'], x['sample'], x['run_label'])):
            writer.writerow(row)

    print('Done.')
    print(f'Output directory: {args.outdir}')
    print('Files created:')
    print(f'  {event_out}')
    print(f'  {collapsed_out}')
    print(f'  {gene_summary_out}')
    print(f'  {suppressed_out}')


if __name__ == '__main__':
    main()
