#!/usr/bin/env python3

import os
import re
import csv
import glob
import argparse
from collections import defaultdict

# ------------------------------------------------------------
# Manual GFF patcher for OsHV-1 merged_safe GFF3 files
# - applies curated per-sample/per-gene replacements from CSV
# - deletes existing features for each targeted gene and replaces them
#   with Gene/mRNA/CDS or misc_feature according to decision
# - handles duplicate genes (Multi=Dual) by adding one replacement per CSV row
# - also renames two incorrect labels globally:
#     MF509813.1_cds_0027 -> ORF24
#     MF509813.1_cds_0043 -> ORF40
# ------------------------------------------------------------

ORF_REGEX = re.compile(r'\bORF[_\-\s]?[A-Za-z0-9.]+\b', re.IGNORECASE)
ALT_OR_REGEX = re.compile(r'^OR(\d+[A-Za-z0-9.]*)$', re.IGNORECASE)

GLOBAL_RENAMES = {
    'MF509813.1_cds_0027': 'ORF24',
    'MF509813.1_cds_0043': 'ORF40',
}


def normalize_gene_name(name):
    if name is None:
        return None
    name = str(name).strip()
    if not name or name.upper() == 'NA':
        return None
    if name in GLOBAL_RENAMES:
        return GLOBAL_RENAMES[name]
    m_alt = ALT_OR_REGEX.match(name)
    if m_alt:
        return f"ORF{m_alt.group(1).upper()}"
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


def attrs_to_str(attrs):
    preferred = [
        'ID', 'Parent', 'Name', 'gene', 'locus_tag', 'product', 'protein_id',
        'partial', 'start_range', 'end_range', 'incomplete_at_5prime', 'incomplete_at_3prime',
        'note'
    ]
    parts = []
    for k in preferred:
        if k in attrs and attrs[k] not in ('', None):
            parts.append(f'{k}={attrs[k]}')
    for k, v in attrs.items():
        if k not in preferred and v not in ('', None):
            parts.append(f'{k}={v}')
    return ';'.join(parts)


def extract_gene_name(attrs):
    for key in ['gene', 'Name', 'locus_tag', 'product', 'note', 'ID']:
        if key in attrs:
            g = normalize_gene_name(attrs[key])
            if g:
                return g
    return None


def parse_gff(gff_file):
    header = []
    features = []
    with open(gff_file) as fh:
        for line in fh:
            if line.startswith('#'):
                header.append(line.rstrip('\n'))
                continue
            toks = line.rstrip('\n').split('\t')
            if len(toks) != 9:
                continue
            seqid, source, ftype, start, end, score, strand, phase, attrs_str = toks
            features.append({
                'seqid': seqid,
                'source': source,
                'type': ftype,
                'start': int(start),
                'end': int(end),
                'score': score,
                'strand': strand,
                'phase': phase,
                'attrs': parse_attrs(attrs_str),
            })
    return header, features


def feature_sort_key(feat):
    order = {'gene': 0, 'mRNA': 1, 'transcript': 1, 'misc_feature': 1, 'exon': 2, 'CDS': 3, 'five_prime_UTR': 4, 'three_prime_UTR': 5}
    return (feat['seqid'], feat['start'], feat['end'], order.get(feat['type'], 99), feat['strand'])


def write_gff(header, features, out_gff):
    with open(out_gff, 'w') as out:
        if header:
            for h in header:
                out.write(h + '\n')
        else:
            out.write('##gff-version 3\n')
        for feat in sorted(features, key=feature_sort_key):
            out.write('\t'.join([
                feat['seqid'], feat['source'], feat['type'], str(feat['start']), str(feat['end']),
                feat['score'], feat['strand'], feat['phase'], attrs_to_str(feat['attrs'])
            ]) + '\n')


def set_partial_attrs(attrs, strand, start, end, incomplete_5=False, incomplete_3=False):
    if incomplete_5 or incomplete_3:
        attrs['partial'] = 'true'
    if incomplete_5:
        attrs['incomplete_at_5prime'] = 'true'
    if incomplete_3:
        attrs['incomplete_at_3prime'] = 'true'
    if strand == '+':
        if incomplete_5:
            attrs['start_range'] = f'.,{start}'
        if incomplete_3:
            attrs['end_range'] = f'{end},.'
    else:
        if incomplete_5:
            attrs['end_range'] = f'{end},.'
        if incomplete_3:
            attrs['start_range'] = f'.,{start}'


def normalize_decision(decision):
    if decision is None:
        return None
    d = str(decision).strip().replace(' ', '_')
    d = d.replace('__', '_')
    aliases = {
        'Partial_end': 'Partial_end',
        'Partial_start': 'Partial_start',
        'Partial_stop': 'Partial_end',
        'Truncation_move': 'Truncation_move',
        'Truncation_stay': 'Truncation_stay',
        'Truncation': 'Truncation_stay',
        'ORF': 'ORF',
        'MISC': 'MISC',
        'Delete': 'Delete',
    }
    # preserve exact aliases regardless of case
    for k, v in aliases.items():
        if d.lower() == k.lower():
            return v
    return d


def parse_orientation(value, fallback=None):
    if value is None:
        return fallback
    v = str(value).strip().lower()
    if v in ('plus', '+'):
        return '+'
    if v in ('minus', '-'):
        return '-'
    return fallback


def used_ids_from_features(features):
    ids = set()
    for f in features:
        if 'ID' in f['attrs']:
            ids.add(f['attrs']['ID'])
    return ids


def make_unique_id(base, used_ids):
    candidate = re.sub(r'[^A-Za-z0-9_.-]+', '_', base)
    if candidate not in used_ids:
        used_ids.add(candidate)
        return candidate
    i = 2
    while True:
        c2 = f'{candidate}_{i}'
        if c2 not in used_ids:
            used_ids.add(c2)
            return c2
        i += 1


def model_gene_bounds(features_for_gene):
    if not features_for_gene:
        return None
    seqid = features_for_gene[0]['seqid']
    strand = features_for_gene[0]['strand']
    start = min(f['start'] for f in features_for_gene)
    end = max(f['end'] for f in features_for_gene)
    return seqid, start, end, strand


def rename_global_labels(features):
    for feat in features:
        attrs = feat['attrs']
        for k, v in list(attrs.items()):
            if v in GLOBAL_RENAMES:
                attrs[k] = GLOBAL_RENAMES[v]
    return features


def collect_features_by_gene(features):
    by_gene = defaultdict(list)
    for feat in features:
        g = extract_gene_name(feat['attrs'])
        if g:
            by_gene[g].append(feat)
    return by_gene


def remove_gene_features(features, target_gene):
    kept = []
    removed = []
    for feat in features:
        g = extract_gene_name(feat['attrs'])
        if g == target_gene:
            removed.append(feat)
        else:
            kept.append(feat)
    return kept, removed


def infer_seqid_from_features(sample, features, removed_features):
    # prefer removed features for the same gene; else the dominant seqid in file; else sample
    seqids = [f['seqid'] for f in removed_features if f.get('seqid')]
    if seqids:
        return seqids[0]
    counts = defaultdict(int)
    for f in features:
        counts[f['seqid']] += 1
    if counts:
        return sorted(counts.items(), key=lambda x: x[1], reverse=True)[0][0]
    return sample


def infer_coords_for_stay(removed_features):
    if not removed_features:
        return None
    b = model_gene_bounds(removed_features)
    if not b:
        return None
    seqid, start, end, strand = b
    # return transcript-ish coordinates: for minus orientation, report larger first for internal use if desired
    if strand == '-':
        return seqid, end, start, strand
    return seqid, start, end, strand


def build_replacement_features(sample, gene, seqid, start_in, end_in, strand, decision, used_ids, row_idx=1):
    # GFF coordinates always ascending
    start = min(start_in, end_in)
    end = max(start_in, end_in)
    base_tag = f'{sample}__{gene}__manual__{decision}__{row_idx}'
    note_base = f'manual_curation_decision={decision}'

    if decision == 'MISC':
        misc_id = make_unique_id(base_tag + '__misc', used_ids)
        attrs = {
            'ID': misc_id,
            'Name': gene,
            'gene': gene,
            'note': f'gene "{gene}" likely present but assembly gaps make unclear;{note_base}'
        }
        return [{
            'seqid': seqid, 'source': 'ManualRedo', 'type': 'misc_feature',
            'start': start, 'end': end, 'score': '.', 'strand': strand, 'phase': '.', 'attrs': attrs
        }]

    gene_id = make_unique_id(base_tag + '__gene', used_ids)
    mrna_id = make_unique_id(base_tag + '__mrna', used_ids)
    cds_id = make_unique_id(base_tag + '__cds1', used_ids)

    product = gene
    gene_attrs = {
        'ID': gene_id,
        'Name': gene,
        'gene': gene,
        'locus_tag': gene,
        'product': product,
        'note': note_base,
    }
    mrna_attrs = {
        'ID': mrna_id,
        'Parent': gene_id,
        'Name': gene,
        'gene': gene,
        'locus_tag': gene,
        'product': product,
        'note': note_base,
    }
    cds_attrs = {
        'ID': cds_id,
        'Parent': mrna_id,
        'Name': gene,
        'gene': gene,
        'locus_tag': gene,
        'product': product,
        'protein_id': gene,
        'note': note_base,
    }

    incomplete_5 = False
    incomplete_3 = False
    if decision == 'Partial_start':
        incomplete_5 = True
        extra = 'missing_start_codon_due_to_assembly_gap'
        gene_attrs['note'] += ';' + extra
        mrna_attrs['note'] += ';' + extra
        cds_attrs['note'] += ';' + extra
    elif decision == 'Partial_end':
        incomplete_3 = True
        extra = 'missing_end_codon_due_to_assembly_gap'
        gene_attrs['note'] += ';' + extra
        mrna_attrs['note'] += ';' + extra
        cds_attrs['note'] += ';' + extra
    elif decision in ('Truncation_move', 'Truncation_stay'):
        extra = 'curated_truncated_feature'
        gene_attrs['note'] += ';' + extra
        mrna_attrs['note'] += ';' + extra
        cds_attrs['note'] += ';' + extra
    elif decision == 'ORF':
        gene_attrs['note'] += ';curated_complete_feature'
        mrna_attrs['note'] += ';curated_complete_feature'
        cds_attrs['note'] += ';curated_complete_feature'

    if incomplete_5 or incomplete_3:
        set_partial_attrs(gene_attrs, strand, start, end, incomplete_5, incomplete_3)
        set_partial_attrs(mrna_attrs, strand, start, end, incomplete_5, incomplete_3)
        set_partial_attrs(cds_attrs, strand, start, end, incomplete_5, incomplete_3)

    gene_feat = {'seqid': seqid, 'source': 'ManualRedo', 'type': 'gene', 'start': start, 'end': end, 'score': '.', 'strand': strand, 'phase': '.', 'attrs': gene_attrs}
    mrna_feat = {'seqid': seqid, 'source': 'ManualRedo', 'type': 'mRNA', 'start': start, 'end': end, 'score': '.', 'strand': strand, 'phase': '.', 'attrs': mrna_attrs}
    cds_feat = {'seqid': seqid, 'source': 'ManualRedo', 'type': 'CDS', 'start': start, 'end': end, 'score': '.', 'strand': strand, 'phase': '0', 'attrs': cds_attrs}
    return [gene_feat, mrna_feat, cds_feat]


def read_csv_rows(csv_path):
    rows = []
    with open(csv_path, newline='') as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            sample = row.get('Sample', '').strip()
            gene = normalize_gene_name(row.get('Gene'))
            if not sample or not gene:
                continue
            decision = normalize_decision(row.get('decision'))
            start_raw = str(row.get('Start', '')).strip()
            end_raw = str(row.get('End', '')).strip()
            start = None if start_raw in ('', 'NA', 'nan', 'None') else int(float(start_raw))
            end = None if end_raw in ('', 'NA', 'nan', 'None') else int(float(end_raw))
            orient = parse_orientation(row.get('Orientation'))
            multi = str(row.get('Multi', '')).strip()
            rows.append({
                'Sample': sample,
                'Gene': gene,
                'Start': start,
                'End': end,
                'Orientation': orient,
                'decision': decision,
                'Multi': multi,
            })
    return rows


def patch_one_file(gff_path, rows_for_sample, out_path):
    header, features = parse_gff(gff_path)
    features = rename_global_labels(features)
    original_features = list(features)
    used_ids = used_ids_from_features(features)

    # group existing features by gene after global renames
    by_gene_initial = collect_features_by_gene(features)

    summary_rows = []

    # process gene groups (for Multi=Dual etc.)
    rows_by_gene = defaultdict(list)
    for row in rows_for_sample:
        rows_by_gene[row['Gene']].append(row)

    for gene, gene_rows in rows_by_gene.items():
        # remove all current features for this gene first
        features, removed_features = remove_gene_features(features, gene)
        existing_info = infer_coords_for_stay(removed_features)
        seqid_fallback = infer_seqid_from_features(rows_for_sample[0]['Sample'], original_features, removed_features)

        row_counter = 0
        for row in gene_rows:
            row_counter += 1
            decision = row['decision']
            orientation = row['Orientation']
            start = row['Start']
            end = row['End']

            if orientation is None and existing_info is not None:
                orientation = existing_info[3]
            if orientation is None:
                orientation = '+'

            if decision == 'Delete':
                summary_rows.append({
                    'sample': row['Sample'], 'gene': gene, 'decision': decision,
                    'status': 'deleted', 'details': 'current feature removed and not replaced'
                })
                continue

            if decision in ('Truncation_stay',) and (start is None or end is None):
                if existing_info is None:
                    summary_rows.append({
                        'sample': row['Sample'], 'gene': gene, 'decision': decision,
                        'status': 'warning_skipped', 'details': 'no existing coordinates available for stay decision'
                    })
                    continue
                _, s_exist, e_exist, strand_exist = existing_info
                orientation = strand_exist
                # existing_info returns transcript-order-esque coords for minus already
                start = s_exist
                end = e_exist
            elif decision in ('Truncation_move', 'Partial_start', 'Partial_end', 'ORF', 'MISC'):
                if start is None or end is None:
                    summary_rows.append({
                        'sample': row['Sample'], 'gene': gene, 'decision': decision,
                        'status': 'warning_skipped', 'details': 'coordinates missing for decision requiring coordinates'
                    })
                    continue
            elif decision not in ('Delete', 'Truncation_stay') and (start is None or end is None):
                # generic fallback for Truncation alias or unknown similar entries
                if existing_info is not None and decision.startswith('Truncation'):
                    _, s_exist, e_exist, strand_exist = existing_info
                    orientation = strand_exist
                    start = s_exist
                    end = e_exist
                else:
                    summary_rows.append({
                        'sample': row['Sample'], 'gene': gene, 'decision': decision,
                        'status': 'warning_skipped', 'details': 'could not determine coordinates'
                    })
                    continue

            seqid = seqid_fallback
            if removed_features:
                seqid = removed_features[0]['seqid']
            replacement_feats = build_replacement_features(
                row['Sample'], gene, seqid, start, end, orientation, decision, used_ids, row_idx=row_counter
            )
            features.extend(replacement_feats)
            summary_rows.append({
                'sample': row['Sample'], 'gene': gene, 'decision': decision,
                'status': 'replaced', 'details': f'added_{len(replacement_feats)}_feature(s)'
            })

    # Also apply global label rename to the remaining non-targeted features already done
    write_gff(header, features, out_path)
    return summary_rows


def main():
    ap = argparse.ArgumentParser(description='Apply manually curated coordinate/decision updates to each .merged_safe.gff3 file using a CSV table.')
    ap.add_argument('--merged-dir', required=True, help='Directory containing per-sample merged_safe.gff3 files (e.g. oshv1_annotation_merged_MF509813plus)')
    ap.add_argument('--csv', required=True, help='CSV file with manual redo instructions')
    ap.add_argument('--outdir', required=True, help='Output directory for patched GFF3 files')
    ap.add_argument('--inplace', action='store_true', help='Overwrite original files under --merged-dir instead of writing to --outdir')
    args = ap.parse_args()

    rows = read_csv_rows(args.csv)
    if not rows:
        raise RuntimeError('No usable rows found in CSV')

    rows_by_sample = defaultdict(list)
    for row in rows:
        rows_by_sample[row['Sample']].append(row)

    merged_patterns = [
        os.path.join(args.merged_dir, '*', '*.merged_safe.gff3'),
        os.path.join(args.merged_dir, '*', '*.gff3'),
        os.path.join(args.merged_dir, '*.gff3'),
    ]
    gff_files = []
    seen = set()
    for pat in merged_patterns:
        for fp in sorted(glob.glob(pat)):
            if fp not in seen:
                gff_files.append(fp)
                seen.add(fp)
    if not gff_files:
        raise RuntimeError(f'No GFF3 files found under {args.merged_dir}')

    os.makedirs(args.outdir, exist_ok=True)
    os.makedirs(os.path.join(args.outdir, 'logs'), exist_ok=True)

    overall_summary = []
    unmatched_samples = set(rows_by_sample.keys())

    for gff in gff_files:
        base = os.path.basename(gff)
        sample = None
        # try parent dir name first, then basename without suffix
        parent = os.path.basename(os.path.dirname(gff))
        stem = base.replace('.merged_safe.gff3', '').replace('.gff3', '')
        if parent in rows_by_sample:
            sample = parent
        elif stem in rows_by_sample:
            sample = stem
        else:
            # still apply global renames to every file
            sample = parent if parent else stem

        out_path = gff if args.inplace else os.path.join(args.outdir, parent if parent else sample, base)
        if not args.inplace:
            os.makedirs(os.path.dirname(out_path), exist_ok=True)

        if sample in rows_by_sample:
            unmatched_samples.discard(sample)
            summary_rows = patch_one_file(gff, rows_by_sample[sample], out_path)
            overall_summary.extend(summary_rows)
        else:
            # no manual rows for this sample; still rename the two incorrect names globally
            header, features = parse_gff(gff)
            features = rename_global_labels(features)
            write_gff(header, features, out_path)
            overall_summary.append({
                'sample': sample, 'gene': '', 'decision': '', 'status': 'copied_with_global_renames_only', 'details': os.path.basename(gff)
            })

    summary_tsv = os.path.join(args.outdir, 'logs', 'manual_redo_summary.tsv')
    with open(summary_tsv, 'w', newline='') as out:
        writer = csv.DictWriter(out, fieldnames=['sample', 'gene', 'decision', 'status', 'details'], delimiter='\t')
        writer.writeheader()
        for row in overall_summary:
            writer.writerow(row)

    unmatched_tsv = os.path.join(args.outdir, 'logs', 'manual_redo_unmatched_samples.tsv')
    with open(unmatched_tsv, 'w', newline='') as out:
        writer = csv.writer(out, delimiter='\t')
        writer.writerow(['sample'])
        for s in sorted(unmatched_samples):
            writer.writerow([s])

    print('Done.')
    print(f'Patched GFF3 files written to: {args.merged_dir if args.inplace else args.outdir}')
    print(f'Summary: {summary_tsv}')
    print(f'Unmatched CSV samples: {unmatched_tsv}')


if __name__ == '__main__':
    main()
