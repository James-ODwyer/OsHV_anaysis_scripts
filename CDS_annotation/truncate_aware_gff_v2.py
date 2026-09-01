#!/usr/bin/env python3

import os
import re
import csv
import glob
import argparse
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq

ORF_REGEX = re.compile(r'\bORF[_\-\s]?[A-Za-z0-9.]+\b', re.IGNORECASE)


def parse_attrs(attr_str):
    attrs = {}
    for item in attr_str.strip().split(';'):
        if '=' in item:
            k, v = item.split('=', 1)
            attrs[k] = v
    return attrs


def attrs_to_str(attrs):
    preferred = ['ID', 'Parent', 'Name', 'gene', 'locus_tag', 'product', 'protein_id', 'partial', 'start_range', 'end_range', 'note']
    pieces = []
    for k in preferred:
        if k in attrs and attrs[k] not in ('', None):
            pieces.append(f'{k}={attrs[k]}')
    for k, v in attrs.items():
        if k not in preferred and v not in ('', None):
            pieces.append(f'{k}={v}')
    return ';'.join(pieces)


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


def extract_gene_name(attrs):
    for key in ['gene', 'Name', 'locus_tag', 'product', 'note', 'ID']:
        if key in attrs:
            val = normalize_gene_name(attrs[key])
            if val:
                return val
    return None


def parse_gff(gff_file):
    header = []
    features = []
    with open(gff_file) as fh:
        for line in fh:
            if line.startswith('#'):
                header.append(line.rstrip('\n'))
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) != 9:
                continue
            seqid, source, ftype, start, end, score, strand, phase, attrs_str = parts
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


def find_models(features):
    gene_models = {}
    mrna_to_gene = {}
    orphans = []

    for feat in features:
        if feat['type'] == 'gene':
            gid = feat['attrs'].get('ID')
            if not gid:
                gid = f"gene__{feat['seqid']}__{feat['start']}__{feat['end']}"
                feat['attrs']['ID'] = gid
            gene_models[gid] = {
                'gene_id': gid,
                'gene': feat,
                'mrnas': [],
                'cds': [],
                'other': []
            }

    for feat in features:
        if feat['type'] in {'mRNA', 'transcript'}:
            parent = feat['attrs'].get('Parent')
            if parent in gene_models:
                gene_models[parent]['mrnas'].append(feat)
                if 'ID' in feat['attrs']:
                    mrna_to_gene[feat['attrs']['ID']] = parent
            else:
                gid = f"orphan_gene_from_{feat['attrs'].get('ID', feat['seqid'])}"
                if gid not in gene_models:
                    gene_models[gid] = {'gene_id': gid, 'gene': None, 'mrnas': [], 'cds': [], 'other': []}
                gene_models[gid]['mrnas'].append(feat)
                if 'ID' in feat['attrs']:
                    mrna_to_gene[feat['attrs']['ID']] = gid

    for feat in features:
        if feat['type'] in {'gene', 'mRNA', 'transcript'}:
            continue
        parent = feat['attrs'].get('Parent')
        gid = None
        if parent in gene_models:
            gid = parent
        elif parent in mrna_to_gene:
            gid = mrna_to_gene[parent]
        if gid is None:
            orphans.append(feat)
            continue
        if feat['type'] == 'CDS':
            gene_models[gid]['cds'].append(feat)
        else:
            gene_models[gid]['other'].append(feat)
    return gene_models, orphans


def transcript_order(cds_list, strand):
    if strand == '+':
        return sorted(cds_list, key=lambda x: (x['start'], x['end']))
    return sorted(cds_list, key=lambda x: (-x['end'], -x['start']))


def feature_nt_seq(genome_seq, feat):
    subseq = genome_seq[feat['seqid']][feat['start'] - 1:feat['end']]
    if feat['strand'] == '-':
        subseq = str(Seq(subseq).reverse_complement())
    return subseq.upper()


def concat_cds_seq(genome_seq, cds_feats):
    if not cds_feats:
        return ''
    ordered = transcript_order(cds_feats, cds_feats[0]['strand'])
    return ''.join(feature_nt_seq(genome_seq, f) for f in ordered)


def leading_n_count(seq):
    i = 0
    while i < len(seq) and seq[i] == 'N':
        i += 1
    return i


def trailing_n_count(seq):
    i = 0
    j = len(seq) - 1
    while j >= 0 and seq[j] == 'N':
        i += 1
        j -= 1
    return i


def trim_end_from_feature(feat, consume, from_5prime):
    length = feat['end'] - feat['start'] + 1
    if consume <= 0:
        return None, 0
    if consume >= length:
        return None, length
    new_feat = dict(feat)
    new_feat['attrs'] = dict(feat['attrs'])
    if feat['strand'] == '+':
        if from_5prime:
            new_feat['start'] = feat['start'] + consume
        else:
            new_feat['end'] = feat['end'] - consume
    else:
        if from_5prime:
            new_feat['end'] = feat['end'] - consume
        else:
            new_feat['start'] = feat['start'] + consume
    return new_feat, consume


def trim_cds_features(cds_feats, trim5=0, trim3=0):
    if not cds_feats:
        return [], 0, 0
    strand = cds_feats[0]['strand']
    ordered = transcript_order(cds_feats, strand)

    remaining = []
    consume = trim5
    for feat in ordered:
        if consume <= 0:
            remaining.append(dict(feat, attrs=dict(feat['attrs'])))
            continue
        new_feat, consumed = trim_end_from_feature(feat, consume, from_5prime=True)
        consume -= consumed
        if new_feat is not None:
            remaining.append(new_feat)
    actual_trim5 = trim5 - max(consume, 0)

    if not remaining:
        return [], actual_trim5, 0

    ordered2 = transcript_order(remaining, strand)
    consume = trim3
    kept_rev = []
    for feat in reversed(ordered2):
        if consume <= 0:
            kept_rev.append(dict(feat, attrs=dict(feat['attrs'])))
            continue
        new_feat, consumed = trim_end_from_feature(feat, consume, from_5prime=False)
        consume -= consumed
        if new_feat is not None:
            kept_rev.append(new_feat)
    actual_trim3 = trim3 - max(consume, 0)
    final_feats = list(reversed(kept_rev))
    return transcript_order(final_feats, strand), actual_trim5, actual_trim3


def set_partial_attrs(feat, start_partial=False, end_partial=False):
    attrs = feat['attrs']
    if start_partial or end_partial:
        attrs['partial'] = 'true'
    if feat['strand'] == '+':
        if start_partial:
            attrs['start_range'] = f'.,{feat["start"]}'
        if end_partial:
            attrs['end_range'] = f'{feat["end"]},.'
    else:
        if start_partial:
            attrs['end_range'] = f'{feat["end"]},.'
        if end_partial:
            attrs['start_range'] = f'.,{feat["start"]}'


def append_note(attrs, text):
    if not text:
        return
    old = attrs.get('note', '')
    attrs['note'] = f'{old};{text}'.strip(';')


def update_parent_bounds(model):
    cds = model['cds']
    if not cds:
        return
    new_start = min(f['start'] for f in cds)
    new_end = max(f['end'] for f in cds)
    if model['gene'] is not None:
        model['gene']['start'] = new_start
        model['gene']['end'] = new_end
    for mrna in model['mrnas']:
        mrna['start'] = new_start
        mrna['end'] = new_end


def first_internal_stop_nt_len(seq):
    usable_len = (len(seq) // 3) * 3
    if usable_len < 3:
        return None
    cds_seq = seq[:usable_len]
    pep = str(Seq(cds_seq).translate(table=1, to_stop=False))
    stop_positions = [i for i, aa in enumerate(pep) if aa == '*']
    if not stop_positions:
        return None
    idx = stop_positions[0]
    if idx == len(pep) - 1:
        return None
    return (idx + 1) * 3


def compute_ref_lengths(reference_gff):
    gene_id_to_name = {}
    mrna_to_gene = {}
    cds_len_by_gene = defaultdict(int)
    with open(reference_gff) as fh:
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
                gid = attrs.get('ID')
                gname = extract_gene_name(attrs)
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
                    cds_len_by_gene[gname] += (end - start + 1)
    return dict(cds_len_by_gene)


def model_gene_name(model):
    if model['gene'] is not None:
        g = extract_gene_name(model['gene']['attrs'])
        if g:
            return g
    for feat in model['mrnas']:
        g = extract_gene_name(feat['attrs'])
        if g:
            return g
    for feat in model['cds']:
        g = extract_gene_name(feat['attrs'])
        if g:
            return g
    return None


def model_gene_id(model):
    if model['gene'] is not None:
        return model['gene']['attrs'].get('ID', model['gene_id'])
    if model['mrnas']:
        return model['mrnas'][0]['attrs'].get('Parent') or model['mrnas'][0]['attrs'].get('ID') or model['gene_id']
    return model['gene_id']


def has_gap_near_start(seq, within_bp=3):
    if not seq:
        return False
    # detect an N immediately after, overlapping, or within 3 bp of the first codon region
    # effectively checks the first 6 bp of transcript sequence
    window_end = min(len(seq), 3 + within_bp)
    return 'N' in seq[:window_end]


def remove_bad_start_gap_model(genome_seq, model, ref_lengths, gap_within_bp=3, min_fraction=0.9):
    gene_name = model_gene_name(model)
    if not gene_name or gene_name not in ref_lengths:
        return False, ''
    if not model['cds']:
        return False, ''
    seq = concat_cds_seq(genome_seq, model['cds'])
    if not seq:
        return False, ''
    observed_nt = sum((f['end'] - f['start'] + 1) for f in model['cds'])
    expected_nt = ref_lengths[gene_name]
    if expected_nt <= 0:
        return False, ''
    small = observed_nt < (expected_nt * min_fraction)
    gap_near = has_gap_near_start(seq, within_bp=gap_within_bp)
    if gap_near and small:
        return True, f'gap_within_{gap_within_bp}bp_of_start_and_shorter_than_{min_fraction:.2f}x_reference;observed_nt={observed_nt};expected_nt={expected_nt}'
    return False, ''


def fix_model(genome_seq, model):
    cds = model['cds']
    gene_name = model_gene_name(model)
    gene_id = model_gene_id(model)
    if not cds:
        return {'gene': gene_name, 'gene_id': gene_id, 'status': 'no_cds', 'trim5_bp': 0, 'trim3_bp': 0, 'internal_stop_trim_bp': 0, 'five_partial': False, 'three_partial': False}

    seq = concat_cds_seq(genome_seq, cds)
    if not seq:
        return {'gene': gene_name, 'gene_id': gene_id, 'status': 'empty_seq', 'trim5_bp': 0, 'trim3_bp': 0, 'internal_stop_trim_bp': 0, 'five_partial': False, 'three_partial': False}

    lead_ns = leading_n_count(seq)
    trail_ns = trailing_n_count(seq)
    notes = []

    trimmed_cds, actual_trim5, actual_trim3 = trim_cds_features(cds, trim5=lead_ns, trim3=trail_ns)
    if actual_trim5 > 0:
        notes.append(f'trimmed_{actual_trim5}bp_from_5prime_due_to_terminal_Ns')
    if actual_trim3 > 0:
        notes.append(f'trimmed_{actual_trim3}bp_from_3prime_due_to_terminal_Ns')
    if not trimmed_cds:
        model['cds'] = []
        return {'gene': gene_name, 'gene_id': gene_id, 'status': 'fully_trimmed_to_zero', 'trim5_bp': actual_trim5, 'trim3_bp': actual_trim3, 'internal_stop_trim_bp': 0, 'five_partial': actual_trim5 > 0, 'three_partial': actual_trim3 > 0}

    seq2 = concat_cds_seq(genome_seq, trimmed_cds)
    stop_nt_len = first_internal_stop_nt_len(seq2)
    actual_stop_trim = 0
    if stop_nt_len is not None and stop_nt_len < len(seq2):
        actual_stop_trim = len(seq2) - stop_nt_len
        trimmed_cds, _, actual_stop_trim = trim_cds_features(trimmed_cds, trim5=0, trim3=actual_stop_trim)
        notes.append(f'truncated_at_first_internal_stop;trimmed_{actual_stop_trim}bp_from_3prime')

    model['cds'] = trimmed_cds
    update_parent_bounds(model)

    five_partial = actual_trim5 > 0
    three_partial = actual_trim3 > 0
    if model['gene'] is not None:
        set_partial_attrs(model['gene'], start_partial=five_partial, end_partial=three_partial)
        for note in notes:
            append_note(model['gene']['attrs'], note)
    for mrna in model['mrnas']:
        set_partial_attrs(mrna, start_partial=five_partial, end_partial=three_partial)
        for note in notes:
            append_note(mrna['attrs'], note)
    for cds_feat in model['cds']:
        set_partial_attrs(cds_feat, start_partial=five_partial, end_partial=three_partial)
        for note in notes:
            append_note(cds_feat['attrs'], note)
        if five_partial:
            cds_feat['phase'] = str((3 - (actual_trim5 % 3)) % 3)

    return {
        'gene': gene_name,
        'gene_id': gene_id,
        'status': 'updated',
        'trim5_bp': actual_trim5,
        'trim3_bp': actual_trim3,
        'internal_stop_trim_bp': actual_stop_trim,
        'five_partial': five_partial,
        'three_partial': three_partial,
    }


def flatten_models(models, removed_ids, orphans):
    feats = []
    for gid, model in models.items():
        if gid in removed_ids:
            continue
        if model['gene'] is not None:
            feats.append(model['gene'])
        feats.extend(model['mrnas'])
        feats.extend(model['other'])
        feats.extend(model['cds'])
    feats.extend(orphans)
    return feats


def feature_sort_key(feat):
    order = {'gene': 0, 'mRNA': 1, 'transcript': 1, 'exon': 2, 'CDS': 3, 'five_prime_UTR': 4, 'three_prime_UTR': 5}
    return (feat['seqid'], feat['start'], feat['end'], order.get(feat['type'], 99), feat['strand'])


def write_gff(header, feats, out_gff):
    with open(out_gff, 'w') as out:
        if header:
            for h in header:
                out.write(h + '\n')
        else:
            out.write('##gff-version 3\n')
        for feat in sorted(feats, key=feature_sort_key):
            out.write('\t'.join([
                feat['seqid'], feat['source'], feat['type'], str(feat['start']), str(feat['end']),
                feat['score'], feat['strand'], feat['phase'], attrs_to_str(feat['attrs'])
            ]) + '\n')


def main():
    ap = argparse.ArgumentParser(description='Make a GFF truncate-aware by trimming CDSs at terminal Ns and first internal stop codon, and optionally remove likely artefactual annotations caused by a large N gap near the gene start.')
    ap.add_argument('--genome', required=True, help='Genome FASTA for the sample corresponding to the GFF')
    ap.add_argument('--gff', required=True, help='Input GFF3')
    ap.add_argument('--out-gff', required=True, help='Output corrected GFF3')
    ap.add_argument('--summary', required=False, help='Optional per-gene summary TSV')
    ap.add_argument('--reference-gff', required=False, help='Primary reference GFF3 for expected CDS lengths; enables gap-near-start removal')
    ap.add_argument('--removed-log', required=False, help='Optional TSV logging annotations removed because a gap of Ns occurs near the start codon and the observed gene is shorter than expected')
    ap.add_argument('--gap-within-bp', type=int, default=3, help='Look for N gaps within this many bp of the start codon (default: 3)')
    ap.add_argument('--min-fraction-of-reference', type=float, default=0.90, help='Remove start-gap cases only if observed CDS length is below this fraction of expected reference CDS length (default: 0.90)')
    args = ap.parse_args()

    genome_seq = {rec.id: str(rec.seq).upper() for rec in SeqIO.parse(args.genome, 'fasta')}
    if not genome_seq:
        raise RuntimeError(f'No FASTA records found in {args.genome}')

    header, feats = parse_gff(args.gff)
    models, orphans = find_models(feats)
    ref_lengths = compute_ref_lengths(args.reference_gff) if args.reference_gff else {}

    summaries = []
    removed_rows = []
    removed_ids = set()

    for gid, model in models.items():
        res = fix_model(genome_seq, model)
        summaries.append(res)
        if args.reference_gff and model['cds']:
            remove, reason = remove_bad_start_gap_model(
                genome_seq, model, ref_lengths,
                gap_within_bp=args.gap_within_bp,
                min_fraction=args.min_fraction_of_reference
            )
            if remove:
                removed_ids.add(gid)
                removed_rows.append({
                    'sample': next(iter(genome_seq.keys())),
                    'gene': model_gene_name(model) or '',
                    'gene_id': model_gene_id(model) or gid,
                    'reason': reason,
                })

    out_feats = flatten_models(models, removed_ids, orphans)
    write_gff(header, out_feats, args.out_gff)

    if args.summary:
        with open(args.summary, 'w', newline='') as out:
            fieldnames = ['gene', 'gene_id', 'status', 'trim5_bp', 'trim3_bp', 'internal_stop_trim_bp', 'five_partial', 'three_partial']
            writer = csv.DictWriter(out, fieldnames=fieldnames, delimiter='\t')
            writer.writeheader()
            for row in summaries:
                writer.writerow(row)

    if args.removed_log:
        os.makedirs(os.path.dirname(args.removed_log), exist_ok=True)
        with open(args.removed_log, 'w', newline='') as out:
            fieldnames = ['sample', 'gene', 'gene_id', 'reason']
            writer = csv.DictWriter(out, fieldnames=fieldnames, delimiter='\t')
            writer.writeheader()
            for row in removed_rows:
                writer.writerow(row)

    print('Done')
    print(f'Wrote corrected GFF: {args.out_gff}')
    if args.summary:
        print(f'Wrote summary TSV: {args.summary}')
    if args.removed_log:
        print(f'Wrote removed-annotation log: {args.removed_log}')


if __name__ == '__main__':
    main()
