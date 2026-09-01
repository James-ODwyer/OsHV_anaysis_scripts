#!/usr/bin/env python3

import os
import re
import csv
import glob
import argparse
from collections import defaultdict

# ------------------------------------------------------------
# Prepare NCBI 5-column feature tables (.tbl) from per-sample
# FASTA + GFF3 for OsHV-1 genomes.
#
# Behaviour:
#   - adds a source feature spanning the whole sequence
#   - keeps only CDS and misc_feature from GFF3
#   - adds assembly_gap features from runs of N in sequence
#   - emits partial CDS using < / > location syntax
#   - sets codon_start from GFF phase for 5' partial CDS when needed
#   - strips internal workflow notes (manual curation / truncation /
#     recovery bookkeeping) so the .tbl is cleaner for GenBank.
# ------------------------------------------------------------

ORF_REGEX = re.compile(r'\bORF[_\-\s]?[A-Za-z0-9.]+\b', re.IGNORECASE)

# Notes matching any of these patterns are omitted from the submission table.
INTERNAL_NOTE_PATTERNS = [
    re.compile(r'^manual_curation_decision=', re.IGNORECASE),
    re.compile(r'^truncated_at_first_internal_stop', re.IGNORECASE),
    re.compile(r'^trimmed_\d+bp_from_[53]prime_due_to_terminal_Ns$', re.IGNORECASE),
    re.compile(r'^missing_start_codon_due_to_assembly_gap$', re.IGNORECASE),
    re.compile(r'^missing_end_codon_due_to_assembly_gap$', re.IGNORECASE),
    re.compile(r'^curated_truncated_feature$', re.IGNORECASE),
    re.compile(r'^curated_complete_feature$', re.IGNORECASE),
    re.compile(r'^rescued_partial_by_tblastn', re.IGNORECASE),
    re.compile(r'^merged_from_other_reference$', re.IGNORECASE),
    re.compile(r'^inferred_from_homology', re.IGNORECASE),
    re.compile(r'^best_hit=', re.IGNORECASE),
    re.compile(r'^pident=', re.IGNORECASE),
    re.compile(r'^qcov=', re.IGNORECASE),
    re.compile(r'^evalue=', re.IGNORECASE),
]


def read_single_fasta(path):
    header = None
    seq_chunks = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip('\n')
            if not line:
                continue
            if line.startswith('>'):
                if header is not None:
                    raise RuntimeError(f'More than one FASTA record found in {path}; expected one per sample file')
                header = line[1:].strip().split()[0]
            else:
                seq_chunks.append(line.strip())
    if header is None:
        raise RuntimeError(f'No FASTA record found in {path}')
    return header, ''.join(seq_chunks).upper()


def parse_attrs(attr_str):
    attrs = {}
    for item in attr_str.strip().split(';'):
        if '=' in item:
            k, v = item.split('=', 1)
            attrs[k] = v
    return attrs


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
    feats = []
    with open(gff_file) as fh:
        for line in fh:
            if line.startswith('#'):
                header.append(line.rstrip('\n'))
                continue
            toks = line.rstrip('\n').split('\t')
            if len(toks) != 9:
                continue
            seqid, source, ftype, start, end, score, strand, phase, attrs_str = toks
            feats.append({
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
    return header, feats


def find_runs_of_ns(seq, min_len=10):
    runs = []
    i = 0
    n = len(seq)
    while i < n:
        if seq[i] == 'N':
            j = i
            while j < n and seq[j] == 'N':
                j += 1
            run_len = j - i
            if run_len >= min_len:
                runs.append((i + 1, j, run_len))
            i = j
        else:
            i += 1
    return runs


def year_from_sample(sample):
    m = re.search(r'(19|20)\d{2}', sample)
    return m.group(0) if m else None


def build_source_feature(seqid, seq_len, sample, organism, mol_type, taxon, host):
    quals = [
        ('organism', organism),
        ('mol_type', mol_type),
        ('db_xref', f'taxon:{taxon}'),
        ('host', host),
        ('isolate', sample),
    ]
    year = year_from_sample(sample)
    if year:
        quals.append(('collection_date', year))
    return {'seqid': seqid, 'start': 1, 'end': seq_len, 'type': 'source', 'quals': quals, 'strand': '+'}


def build_gap_features(seqid, seq, gap_min, gap_type, linkage_evidence):
    feats = []
    for s, e, length in find_runs_of_ns(seq, min_len=gap_min):
        feats.append({
            'seqid': seqid,
            'start': s,
            'end': e,
            'type': 'assembly_gap',
            'quals': [
                ('estimated_length', str(length)),
                ('gap_type', gap_type),
                ('linkage_evidence', linkage_evidence),
            ],
            'strand': '+'
        })
    return feats


def is_partial_5(attrs):
    return attrs.get('incomplete_at_5prime', '').lower() == 'true'


def is_partial_3(attrs):
    return attrs.get('incomplete_at_3prime', '').lower() == 'true'


def location_strings_for_tbl(start, end, strand, partial5=False, partial3=False):
    # For the 5-column feature table, minus-strand features are represented by
    # reversing the coordinates. Partial symbols: '<' always column 1, '>' always column 2.
    if strand == '-':
        col1 = str(end)
        col2 = str(start)
        if partial5:
            col1 = '<' + col1
        if partial3:
            col2 = '>' + col2
    else:
        col1 = str(start)
        col2 = str(end)
        if partial5:
            col1 = '<' + col1
        if partial3:
            col2 = '>' + col2
    return col1, col2


def codon_start_from_phase(phase):
    if phase in ('.', '', None):
        return 1
    try:
        p = int(phase)
    except ValueError:
        return 1
    return p + 1 if p in (0, 1, 2) else 1


def is_internal_note(note_text):
    note_text = str(note_text).strip()
    if not note_text:
        return True
    for pattern in INTERNAL_NOTE_PATTERNS:
        if pattern.search(note_text):
            return True
    return False


def cleaned_notes(raw_note):
    if raw_note in ('', None):
        return []
    notes = []
    for note in str(raw_note).split(';'):
        note = note.strip()
        if not note:
            continue
        if is_internal_note(note):
            continue
        notes.append(note)
    return notes


def convert_gff_feats_to_tbl(seqid, feats, keep_cds_notes=False):
    out = []
    for f in feats:
        if f['type'] not in {'CDS', 'misc_feature'}:
            continue
        attrs = f['attrs']
        if f['type'] == 'misc_feature':
            col1, col2 = location_strings_for_tbl(f['start'], f['end'], f['strand'], False, False)
            quals = []
            note_list = cleaned_notes(attrs.get('note', ''))
            for note in note_list:
                quals.append(('note', note))
            if 'standard_name' in attrs and attrs['standard_name']:
                quals.append(('standard_name', attrs['standard_name']))
            # if no note was preserved and the misc_feature is clearly a manual uncertain region,
            # fall back to a simple gene-based note (only if product/gene info exists and not internal).
            if not note_list:
                g = extract_gene_name(attrs)
                if g and g not in ('None', ''):
                    # only add this fallback for misc_feature, never for CDS.
                    pass
            out.append({'seqid': seqid, 'start': col1, 'end': col2, 'type': 'misc_feature', 'quals': quals, 'strand': f['strand']})
        else:  # CDS
            partial5 = is_partial_5(attrs)
            partial3 = is_partial_3(attrs)
            col1, col2 = location_strings_for_tbl(f['start'], f['end'], f['strand'], partial5, partial3)
            quals = []
            if partial5:
                codon_start = codon_start_from_phase(f['phase'])
                if codon_start != 1:
                    quals.append(('codon_start', str(codon_start)))
            elif 'codon_start' in attrs and attrs['codon_start'] not in ('', None, '1'):
                quals.append(('codon_start', attrs['codon_start']))
            elif f['phase'] not in ('.', '', None):
                cs = codon_start_from_phase(f['phase'])
                if cs != 1:
                    quals.append(('codon_start', str(cs)))

            product = attrs.get('product') or extract_gene_name(attrs) or 'hypothetical protein'
            quals.append(('product', product))
            if 'standard_name' in attrs and attrs['standard_name']:
                quals.append(('standard_name', attrs['standard_name']))
            # By default, omit CDS notes from workflow internals and also suppress CDS notes entirely.
            if keep_cds_notes:
                for note in cleaned_notes(attrs.get('note', '')):
                    quals.append(('note', note))
            out.append({'seqid': seqid, 'start': col1, 'end': col2, 'type': 'CDS', 'quals': quals, 'strand': f['strand']})
    return out


def write_tbl(out_tbl, seqid, features):
    with open(out_tbl, 'w') as out:
        out.write(f'>Features\t{seqid}\n')
        for f in features:
            out.write(f"{f['start']}\t{f['end']}\t{f['type']}\n")
            for qk, qv in f['quals']:
                out.write(f"\t\t\t{qk}\t{qv}\n")


def find_gff_files(input_dir):
    patterns = [
        os.path.join(input_dir, '*', '*.merged_safe.gff3'),
        os.path.join(input_dir, '*', '*.gff3'),
        os.path.join(input_dir, '*.gff3'),
    ]
    found = []
    seen = set()
    for pat in patterns:
        for fp in sorted(glob.glob(pat)):
            if fp not in seen:
                found.append(fp)
                seen.add(fp)
    return found


def guess_sample_name(gff_path):
    base = os.path.basename(gff_path)
    parent = os.path.basename(os.path.dirname(gff_path))
    for suffix in ('.merged_safe.gff3', '.merged_from_multiref.gff3', '.final.gff3', '.gff3'):
        if base.endswith(suffix):
            stem = base[:-len(suffix)]
            if stem:
                return stem
    return parent if parent else os.path.splitext(base)[0]


def main():
    ap = argparse.ArgumentParser(description='Generate NCBI 5-column feature table files (.tbl) from per-sample GFF3 + FASTA, adding source and assembly_gap features and keeping only CDS/misc_feature annotations while stripping internal workflow notes.')
    ap.add_argument('--gff-dir', required=True, help='Directory containing per-sample GFF3 files')
    ap.add_argument('--fasta-dir', required=True, help='Directory containing per-sample FASTA files (sample.fasta)')
    ap.add_argument('--outdir', required=True, help='Output directory for .tbl files')
    ap.add_argument('--gap-min', type=int, default=10, help='Minimum run of Ns to annotate as assembly_gap (default: 10)')
    ap.add_argument('--organism', default='Ostreid herpesvirus 1')
    ap.add_argument('--mol-type', default='genomic DNA')
    ap.add_argument('--taxon', default='26193')
    ap.add_argument('--host', default='Pacific oyster')
    ap.add_argument('--gap-type', default='within scaffold')
    ap.add_argument('--linkage-evidence', default='paired-ends')
    ap.add_argument('--keep-cds-notes', action='store_true', help='If set, keep cleaned non-internal CDS notes. By default CDS notes are omitted.')
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    os.makedirs(os.path.join(args.outdir, 'logs'), exist_ok=True)

    gff_files = find_gff_files(args.gff_dir)
    if not gff_files:
        raise RuntimeError(f'No GFF3 files found under {args.gff_dir}')

    summary_rows = []

    for gff in gff_files:
        sample = guess_sample_name(gff)
        fasta = os.path.join(args.fasta_dir, f'{sample}.fasta')
        if not os.path.isfile(fasta):
            summary_rows.append({'sample': sample, 'status': 'missing_fasta', 'details': fasta})
            continue
        seqid, seq = read_single_fasta(fasta)
        _, gff_feats = parse_gff(gff)
        converted = convert_gff_feats_to_tbl(seqid, gff_feats, keep_cds_notes=args.keep_cds_notes)
        source_feat = build_source_feature(seqid, len(seq), sample, args.organism, args.mol_type, args.taxon, args.host)
        gap_feats = build_gap_features(seqid, seq, args.gap_min, args.gap_type, args.linkage_evidence)

        all_tbl_feats = [source_feat] + gap_feats + converted

        def sort_key_tbl(f):
            def num(x):
                s = str(x).lstrip('<>').strip()
                try:
                    return int(s)
                except ValueError:
                    return 10**12
            order = {'source': 0, 'assembly_gap': 1, 'misc_feature': 2, 'CDS': 3}
            return (num(f['start']), num(f['end']), order.get(f['type'], 99))

        all_tbl_feats = sorted(all_tbl_feats, key=sort_key_tbl)
        out_tbl = os.path.join(args.outdir, f'{sample}.tbl')
        write_tbl(out_tbl, seqid, all_tbl_feats)

        summary_rows.append({
            'sample': sample,
            'status': 'written',
            'details': out_tbl,
            'seq_length': len(seq),
            'assembly_gaps': len(gap_feats),
            'misc_features': sum(1 for x in converted if x['type'] == 'misc_feature'),
            'cds_features': sum(1 for x in converted if x['type'] == 'CDS'),
        })

    summary_tsv = os.path.join(args.outdir, 'logs', 'feature_table_generation_summary.tsv')
    with open(summary_tsv, 'w', newline='') as out:
        fieldnames = ['sample', 'status', 'details', 'seq_length', 'assembly_gaps', 'misc_features', 'cds_features']
        writer = csv.DictWriter(out, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()
        for row in summary_rows:
            writer.writerow(row)

    print('Done.')
    print(f'Feature tables written to: {args.outdir}')
    print(f'Summary: {summary_tsv}')


if __name__ == '__main__':
    main()
