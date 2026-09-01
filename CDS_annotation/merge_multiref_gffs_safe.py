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


def attrs_to_str(attrs):
    preferred = ['ID', 'Parent', 'Name', 'gene', 'locus_tag', 'product', 'protein_id', 'note']
    parts = []
    for k in preferred:
        if k in attrs and attrs[k] not in (None, ''):
            parts.append(f'{k}={attrs[k]}')
    for k, v in attrs.items():
        if k not in preferred and v not in (None, ''):
            parts.append(f'{k}={v}')
    return ';'.join(parts)


def extract_gene_name(attrs):
    for key in ['gene', 'Name', 'locus_tag', 'product', 'note', 'ID']:
        if key in attrs:
            g = normalize_gene_name(attrs[key])
            if g:
                return g
    return None


def find_final_gffs(run_dir):
    pattern = os.path.join(run_dir, '05_final', '*', '*.final.gff3')
    return {os.path.basename(g).replace('.final.gff3', ''): g for g in sorted(glob.glob(pattern))}


def parse_gff(gff_file):
    header, features = [], []
    with open(gff_file) as fh:
        for line in fh:
            if line.startswith('#'):
                header.append(line.rstrip('\n'))
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) != 9:
                continue
            seqid, source, ftype, start, end, score, strand, phase, attr = parts
            features.append({
                'seqid': seqid, 'source': source, 'type': ftype,
                'start': int(start), 'end': int(end), 'score': score,
                'strand': strand, 'phase': phase, 'attrs': parse_attrs(attr)
            })
    return header, features


def build_gene_models(features):
    gene_models = {}
    gene_order = []
    mrna_to_gene = {}

    for feat in features:
        attrs = feat['attrs']
        if feat['type'] == 'gene':
            gid = attrs.get('ID') or f"gene_{feat['seqid']}_{feat['start']}_{feat['end']}"
            if gid not in gene_models:
                gene_models[gid] = {
                    'gene_id': gid,
                    'gene_name': extract_gene_name(attrs),
                    'seqid': feat['seqid'],
                    'start': feat['start'],
                    'end': feat['end'],
                    'strand': feat['strand'],
                    'features': []
                }
                gene_order.append(gid)
            gene_models[gid]['features'].append(feat)
        elif feat['type'] in {'mRNA', 'transcript'}:
            parent = attrs.get('Parent')
            mid = attrs.get('ID')
            if parent and parent in gene_models:
                gene_models[parent]['features'].append(feat)
                if mid:
                    mrna_to_gene[mid] = parent

    for feat in features:
        if feat['type'] not in {'CDS', 'exon', 'five_prime_UTR', 'three_prime_UTR', 'mRNA', 'transcript'}:
            continue
        attrs = feat['attrs']
        assigned_gene = None
        parent = attrs.get('Parent')
        if parent:
            if parent in gene_models:
                assigned_gene = parent
            elif parent in mrna_to_gene:
                assigned_gene = mrna_to_gene[parent]
        if assigned_gene is None:
            gene_name = extract_gene_name(attrs)
            gid = f"orphan_{feat['seqid']}_{feat['start']}_{feat['end']}_{feat['type']}"
            if gid not in gene_models:
                gene_models[gid] = {
                    'gene_id': gid,
                    'gene_name': gene_name,
                    'seqid': feat['seqid'],
                    'start': feat['start'],
                    'end': feat['end'],
                    'strand': feat['strand'],
                    'features': []
                }
                gene_order.append(gid)
            gene_models[gid]['features'].append(feat)
        else:
            if feat not in gene_models[assigned_gene]['features']:
                gene_models[assigned_gene]['features'].append(feat)
            gene_models[assigned_gene]['start'] = min(gene_models[assigned_gene]['start'], feat['start'])
            gene_models[assigned_gene]['end'] = max(gene_models[assigned_gene]['end'], feat['end'])
            if not gene_models[assigned_gene]['gene_name']:
                gene_models[assigned_gene]['gene_name'] = extract_gene_name(attrs)

    by_name = {}
    for gid in gene_order:
        model = gene_models[gid]
        if not model['gene_name']:
            for feat in model['features']:
                g = extract_gene_name(feat['attrs'])
                if g:
                    model['gene_name'] = g
                    break
        if model['gene_name'] and model['gene_name'] not in by_name:
            by_name[model['gene_name']] = model
    return gene_models, gene_order, by_name


def collect_feature_gene_names(features):
    names = set()
    for feat in features:
        if feat['type'] in {'gene', 'mRNA', 'CDS', 'exon'}:
            g = extract_gene_name(feat['attrs'])
            if g:
                names.add(g)
    return names


def collect_intervals(features):
    intervals = defaultdict(list)
    for feat in features:
        if feat['type'] in {'gene', 'mRNA', 'CDS', 'exon'}:
            intervals[(feat['seqid'], feat['strand'])].append((feat['start'], feat['end'], feat['type']))
    return intervals


def overlap_len(a1, a2, b1, b2):
    return max(0, min(a2, b2) - max(a1, b1) + 1)


def model_overlaps_existing(model, base_intervals, min_overlap_fraction=0.10):
    key = (model['seqid'], model['strand'])
    donor_len = model['end'] - model['start'] + 1
    for s, e, ftype in base_intervals.get(key, []):
        ov = overlap_len(model['start'], model['end'], s, e)
        if ov <= 0:
            continue
        frac = ov / donor_len
        if frac >= min_overlap_fraction:
            return True, ov, frac, s, e, ftype
    return False, 0, 0.0, None, None, None


def copy_feature_with_renamed_ids(feat, rename_map):
    new_feat = {
        'seqid': feat['seqid'], 'source': feat['source'], 'type': feat['type'],
        'start': feat['start'], 'end': feat['end'], 'score': feat['score'],
        'strand': feat['strand'], 'phase': feat['phase'], 'attrs': dict(feat['attrs'])
    }
    attrs = new_feat['attrs']
    if 'ID' in attrs and attrs['ID'] in rename_map:
        attrs['ID'] = rename_map[attrs['ID']]
    if 'Parent' in attrs:
        attrs['Parent'] = ','.join(rename_map.get(x, x) for x in attrs['Parent'].split(','))
    attrs['note'] = (attrs.get('note', '') + ';merged_from_other_reference').strip(';')
    return new_feat


def make_unique_id(old_id, sample, gene_name, run_label, used_ids):
    gene_tag = re.sub(r'[^A-Za-z0-9_.-]+', '_', gene_name or 'GENE')
    tail = re.sub(r'[^A-Za-z0-9_.-]+', '_', old_id or 'feature')
    base = f'{sample}__{gene_tag}__from_{run_label}__{tail}'
    candidate = base
    i = 1
    while candidate in used_ids:
        i += 1
        candidate = f'{base}_{i}'
    used_ids.add(candidate)
    return candidate


def feature_sort_key(feat):
    order = {'gene': 0, 'mRNA': 1, 'transcript': 1, 'exon': 2, 'CDS': 3, 'five_prime_UTR': 4, 'three_prime_UTR': 5}
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


def main():
    ap = argparse.ArgumentParser(description='Safely merge missing genes from other reference-based runs without extending truncated genes already present in the base run.')
    ap.add_argument('--base-run', required=True)
    ap.add_argument('--add-run', action='append', required=True, help='LABEL=PATH; may be used multiple times')
    ap.add_argument('--outdir', required=True)
    ap.add_argument('--exclude-genes', default='ORFIN.1,ORFIN.2,ORFIN.3,ORFIN.4')
    ap.add_argument('--min-overlap-fraction', type=float, default=0.10, help='Skip donor gene if it overlaps any existing base feature on the same strand by at least this fraction of donor span (default 0.10)')
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    os.makedirs(os.path.join(args.outdir, 'logs'), exist_ok=True)

    excluded = {normalize_gene_name(x) for x in args.exclude_genes.split(',') if x.strip()}
    base_gffs = find_final_gffs(args.base_run)
    if not base_gffs:
        raise RuntimeError(f'No base final GFF files found under {args.base_run}')

    add_runs = []
    for spec in args.add_run:
        if '=' not in spec:
            raise ValueError(f'Bad --add-run value: {spec}')
        label, path = spec.split('=', 1)
        add_runs.append((label.strip(), path.strip(), find_final_gffs(path.strip())))

    summary_path = os.path.join(args.outdir, 'logs', 'merge_summary_safe.tsv')
    with open(summary_path, 'w', newline='') as summary_fh:
        writer = csv.writer(summary_fh, delimiter='\t')
        writer.writerow(['sample', 'status', 'run_label', 'gene', 'detail'])

        for sample, base_path in sorted(base_gffs.items()):
            header, base_features = parse_gff(base_path)
            _, _, _base_by_name = build_gene_models(base_features)
            base_gene_names = collect_feature_gene_names(base_features)
            base_intervals = collect_intervals(base_features)

            used_ids = {feat['attrs']['ID'] for feat in base_features if 'ID' in feat['attrs']}
            merged_features = list(base_features)

            for run_label, run_path, run_gffs in add_runs:
                donor_path = run_gffs.get(sample)
                if not donor_path:
                    writer.writerow([sample, 'skipped_no_sample', run_label, '', donor_path])
                    continue
                _, donor_features = parse_gff(donor_path)
                _, _, donor_by_name = build_gene_models(donor_features)

                for gene_name, model in sorted(donor_by_name.items()):
                    if not gene_name:
                        continue
                    if gene_name in excluded:
                        writer.writerow([sample, 'skipped_excluded', run_label, gene_name, 'excluded_gene'])
                        continue
                    if gene_name in base_gene_names:
                        writer.writerow([sample, 'skipped_present_in_base', run_label, gene_name, 'gene_name_already_present'])
                        continue

                    overlaps, ov, frac, s, e, ftype = model_overlaps_existing(model, base_intervals, args.min_overlap_fraction)
                    if overlaps:
                        detail = f'overlaps_existing_{ftype}:{s}-{e};ov={ov};frac={frac:.3f}'
                        writer.writerow([sample, 'skipped_overlap', run_label, gene_name, detail])
                        continue

                    old_ids = [feat['attrs']['ID'] for feat in model['features'] if 'ID' in feat['attrs']]
                    rename_map = {old_id: make_unique_id(old_id, sample, gene_name, run_label, used_ids) for old_id in old_ids}
                    copied = [copy_feature_with_renamed_ids(feat, rename_map) for feat in model['features']]
                    merged_features.extend(copied)
                    base_gene_names.add(gene_name)
                    # update intervals to prevent duplicate import from later runs
                    for feat in copied:
                        if feat['type'] in {'gene', 'mRNA', 'CDS', 'exon'}:
                            base_intervals[(feat['seqid'], feat['strand'])].append((feat['start'], feat['end'], feat['type']))
                    writer.writerow([sample, 'added', run_label, gene_name, donor_path])

            out_sample_dir = os.path.join(args.outdir, sample)
            os.makedirs(out_sample_dir, exist_ok=True)
            out_path = os.path.join(out_sample_dir, f'{sample}.merged_safe.gff3')
            write_gff(header, merged_features, out_path)

    print('Done')
    print(f'Output directory: {args.outdir}')
    print(f'Summary log: {summary_path}')


if __name__ == '__main__':
    main()
