#!/usr/bin/env python3
import os
import re
import csv
import glob
import argparse
import subprocess
import tempfile
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
    preferred = [
        'ID', 'Parent', 'Name', 'gene', 'locus_tag', 'product', 'protein_id',
        'partial', 'start_range', 'end_range',
        'incomplete_at_5prime', 'incomplete_at_3prime',
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


def read_fasta_records(path):
    records = []
    header = None
    seq_chunks = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip('\n')
            if not line:
                continue
            if line.startswith('>'):
                if header is not None:
                    records.append((header, ''.join(seq_chunks)))
                header = line[1:].strip()
                seq_chunks = []
            else:
                seq_chunks.append(line.strip())
    if header is not None:
        records.append((header, ''.join(seq_chunks)))
    return records


def parse_reference_proteins(ref_faa):
    by_gene = defaultdict(list)
    for header, seq in read_fasta_records(ref_faa):
        toks = header.split('|')
        accession = toks[0] if len(toks) > 0 else 'NA'
        feature_idx = toks[1] if len(toks) > 1 else 'NA'
        locus_tag = toks[2] if len(toks) > 2 else 'NA'
        gene = toks[3] if len(toks) > 3 else 'NA'
        orf_name = toks[4] if len(toks) > 4 else gene
        protein_id = toks[5] if len(toks) > 5 else 'NA'
        product = toks[6] if len(toks) > 6 else 'hypothetical_protein'
        gene_name = normalize_gene_name(orf_name if orf_name not in ('NA', '') else gene)
        if not gene_name:
            continue
        by_gene[gene_name].append({
            'header': header,
            'seq': seq,
            'accession': accession,
            'feature_idx': feature_idx,
            'locus_tag': locus_tag,
            'gene': gene,
            'orf_name': gene_name,
            'protein_id': protein_id,
            'product': product,
            'qlen': len(seq),
        })
    return by_gene


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


def build_models(features):
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


def find_sample_name(genome_fasta):
    recs = read_fasta_records(genome_fasta)
    if not recs:
        raise RuntimeError(f'No FASTA records found in {genome_fasta}')
    return recs[0][0].split()[0]


def parse_removed_log(removed_log, sample_name=None):
    removed = []
    if not removed_log or not os.path.isfile(removed_log):
        return removed
    with open(removed_log, newline='') as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        for row in reader:
            sample = row.get('sample', '')
            if sample_name and sample and sample != sample_name:
                continue
            removed.append({
                'sample': sample or sample_name or '',
                'gene': normalize_gene_name(row.get('gene')),
                'gene_id': row.get('gene_id', ''),
                'reason': row.get('reason', ''),
            })
    return removed


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


def remove_models_from_gff(features, removed_entries):
    models, orphans = build_models(features)
    removed_gene_names = set(x['gene'] for x in removed_entries if x.get('gene'))
    removed_model_ids = set(x['gene_id'] for x in removed_entries if x.get('gene_id'))
    kept = []
    removed_details = []
    for gid, model in models.items():
        gene_name = model_gene_name(model)
        gene_id = model_gene_id(model)
        if gid in removed_model_ids or gene_id in removed_model_ids or (gene_name and gene_name in removed_gene_names):
            removed_details.append({'gene': gene_name, 'gene_id': gene_id})
            continue
        if model['gene'] is not None:
            kept.append(model['gene'])
        kept.extend(model['mrnas'])
        kept.extend(model['other'])
        kept.extend(model['cds'])
    kept.extend(orphans)
    return kept, removed_details


def make_blast_db(genome_fasta, workdir):
    db_prefix = os.path.join(workdir, 'genome_db')
    cmd = ['makeblastdb', '-in', genome_fasta, '-dbtype', 'nucl', '-out', db_prefix]
    subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    return db_prefix


def write_query_fasta(records, path):
    with open(path, 'w') as out:
        for rec in records:
            out.write(f">{rec['header']}\n")
            seq = rec['seq']
            for i in range(0, len(seq), 80):
                out.write(seq[i:i+80] + '\n')


def run_tblastn(query_faa, db_prefix, out_tsv, threads=1, evalue=1e-5):
    outfmt = '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen'
    cmd = [
        'tblastn', '-query', query_faa, '-db', db_prefix, '-out', out_tsv,
        '-evalue', str(evalue), '-max_target_seqs', '200', '-num_threads', str(threads),
        '-seg', 'no', '-comp_based_stats', '0', '-outfmt', outfmt
    ]
    subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)


def parse_tblastn_hits(tsv_file):
    hits = []
    if not os.path.isfile(tsv_file):
        return hits
    with open(tsv_file) as fh:
        for line in fh:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 13:
                continue
            qseqid, sseqid, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore, qlen = parts[:13]
            qstart = int(qstart); qend = int(qend); sstart = int(sstart); send = int(send); qlen = int(qlen)
            length = int(length); pident = float(pident); bitscore = float(bitscore)
            strand = '+' if sstart <= send else '-'
            start = min(sstart, send);
            end = max(sstart, send)
            hits.append({
                'qseqid': qseqid,
                'sseqid': sseqid,
                'pident': pident,
                'length': length,
                'qstart': qstart,
                'qend': qend,
                'sstart': sstart,
                'send': send,
                'evalue': evalue,
                'bitscore': bitscore,
                'qlen': qlen,
                'strand': strand,
                'start': start,
                'end': end,
            })
    return hits


def merge_intervals(intervals):
    if not intervals:
        return []
    intervals = sorted(intervals)
    merged = [list(intervals[0])]
    for s, e in intervals[1:]:
        if s <= merged[-1][1] + 1:
            merged[-1][1] = max(merged[-1][1], e)
        else:
            merged.append([s, e])
    return [(x[0], x[1]) for x in merged]


def union_len(intervals):
    return sum(e - s + 1 for s, e in merge_intervals(intervals))


def chain_hits_for_group(hits, max_subject_gap=200000, max_query_backtrack=5):
    if not hits:
        return None
    strand = hits[0]['strand']
    if strand == '+':
        hits = sorted(hits, key=lambda h: (h['qstart'], h['start'], -h['bitscore']))
    else:
        hits = sorted(hits, key=lambda h: (h['qstart'], -h['end'], -h['bitscore']))

    chain = []
    last_qend = 0
    last_subject_anchor = None
    for h in hits:
        if chain:
            if h['qstart'] < last_qend - max_query_backtrack:
                continue
            if strand == '+':
                if h['start'] < last_subject_anchor:
                    continue
                if h['start'] - last_subject_anchor > max_subject_gap:
                    continue
            else:
                if h['end'] > last_subject_anchor:
                    continue
                if last_subject_anchor - h['end'] > max_subject_gap:
                    continue
        chain.append(h)
        last_qend = max(last_qend, h['qend'])
        last_subject_anchor = h['end'] if strand == '+' else h['start']

    q_intervals = merge_intervals([(h['qstart'], h['qend']) for h in chain])
    qcov = union_len(q_intervals) / chain[0]['qlen'] if chain and chain[0]['qlen'] else 0.0
    bitscore = sum(h['bitscore'] for h in chain)
    mean_pident = sum(h['pident'] * h['length'] for h in chain) / max(1, sum(h['length'] for h in chain))
    blocks = merge_intervals([(h['start'], h['end']) for h in chain])
    return {
        'hits': chain,
        'blocks': blocks,
        'q_intervals': q_intervals,
        'qcov': qcov,
        'bitscore': bitscore,
        'pident': mean_pident,
        'qstart': min(h['qstart'] for h in chain),
        'qend': max(h['qend'] for h in chain),
        'qlen': chain[0]['qlen'],
        'seqid': chain[0]['sseqid'],
        'strand': chain[0]['strand'],
        'start': min(b[0] for b in blocks),
        'end': max(b[1] for b in blocks),
    }


def choose_best_partial_hit(hits, min_qcov=0.50):
    groups = defaultdict(list)
    for h in hits:
        groups[(h['qseqid'], h['sseqid'], h['strand'])].append(h)
    chains = []
    for key, grp in groups.items():
        chain = chain_hits_for_group(grp)
        if chain is not None:
            chains.append(chain)
    chains = [c for c in chains if c['qcov'] >= min_qcov]
    if not chains:
        return None
    chains.sort(key=lambda c: (c['qcov'], c['bitscore'], c['pident']), reverse=True)
    return chains[0]


def set_partial_boundary_attrs(attrs, strand, start, end, incomplete_5, incomplete_3):
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


def make_unique_ids(sample_name, gene_name, existing_ids):
    base = re.sub(r'[^A-Za-z0-9_.-]+', '_', f'{sample_name}__{gene_name}__partial')
    i = 1
    while True:
        suffix = '' if i == 1 else f'_{i}'
        gene_id = f'{base}__gene{suffix}'
        mrna_id = f'{base}__mrna{suffix}'
        if gene_id not in existing_ids and mrna_id not in existing_ids:
            existing_ids.update([gene_id, mrna_id])
            return gene_id, mrna_id
        i += 1


def build_partial_features(sample_name, gene_name, ref_meta, chain, existing_ids):
    qstart = chain['qstart']
    qend = chain['qend']
    qlen = chain['qlen']
    start = chain['start']
    end = chain['end']
    strand = chain['strand']
    seqid = chain['seqid']
    blocks = sorted(chain['blocks']) if strand == '+' else sorted(chain['blocks'], reverse=True)
    incomplete_5 = qstart > 1
    incomplete_3 = qend < qlen
    gene_id, mrna_id = make_unique_ids(sample_name, gene_name, existing_ids)
    note = (
        f"rescued_partial_by_chained_tblastn;reference_hit={ref_meta['accession']}|{ref_meta['locus_tag']}|{ref_meta['orf_name']};"
        f"qcov={chain['qcov']:.3f};pident={chain['pident']:.2f};qstart={qstart};qend={qend};qlen={qlen};"
        f"blocks={','.join(f'{b[0]}-{b[1]}' for b in blocks)}"
    )

    gene_attrs = {
        'ID': gene_id,
        'Name': gene_name,
        'gene': gene_name,
        'locus_tag': gene_name,
        'product': ref_meta['product'],
        'note': note,
    }
    mrna_attrs = {
        'ID': mrna_id,
        'Parent': gene_id,
        'Name': gene_name,
        'gene': gene_name,
        'locus_tag': gene_name,
        'product': ref_meta['product'],
        'note': note,
    }
    set_partial_boundary_attrs(gene_attrs, strand, start, end, incomplete_5, incomplete_3)
    set_partial_boundary_attrs(mrna_attrs, strand, start, end, incomplete_5, incomplete_3)

    gene_feat = {'seqid': seqid, 'source': 'PartialTblastnChain', 'type': 'gene', 'start': start, 'end': end, 'score': '.', 'strand': strand, 'phase': '.', 'attrs': gene_attrs}
    mrna_feat = {'seqid': seqid, 'source': 'PartialTblastnChain', 'type': 'mRNA', 'start': start, 'end': end, 'score': '.', 'strand': strand, 'phase': '.', 'attrs': mrna_attrs}

    cds_feats = []
    consumed = 0
    block_iter = blocks if strand == '+' else blocks
    for idx, (bstart, bend) in enumerate(block_iter, start=1):
        cds_id = f'{mrna_id}.cds{idx}'
        while cds_id in existing_ids:
            idx += 1
            cds_id = f'{mrna_id}.cds{idx}'
        existing_ids.add(cds_id)
        phase = str(consumed % 3)
        attrs = {
            'ID': cds_id,
            'Parent': mrna_id,
            'Name': gene_name,
            'gene': gene_name,
            'locus_tag': gene_name,
            'product': ref_meta['product'],
            'protein_id': ref_meta['protein_id'],
            'note': note,
        }
        set_partial_boundary_attrs(attrs, strand, bstart, bend, incomplete_5, incomplete_3)
        cds_feats.append({'seqid': seqid, 'source': 'PartialTblastnChain', 'type': 'CDS', 'start': bstart, 'end': bend, 'score': '.', 'strand': strand, 'phase': phase, 'attrs': attrs})
        consumed += (bend - bstart + 1)

    log_row = {
        'gene': gene_name,
        'seqid': seqid,
        'start': start,
        'end': end,
        'strand': strand,
        'incomplete_at_5prime': incomplete_5,
        'incomplete_at_3prime': incomplete_3,
        'qstart': qstart,
        'qend': qend,
        'qlen': qlen,
        'qcov': chain['qcov'],
        'pident': chain['pident'],
        'ref_accession': ref_meta['accession'],
        'ref_locus_tag': ref_meta['locus_tag'],
        'ref_product': ref_meta['product'],
        'blocks': ','.join(f'{b[0]}-{b[1]}' for b in blocks),
    }
    return [gene_feat, mrna_feat] + cds_feats, log_row


def main():
    ap = argparse.ArgumentParser(description='Remove gap-shortened annotations from a GFF and rescue partial gene models by chained tblastn protein homology, allowing split partial CDS features across assembly gaps and annotating incomplete_at_5prime / incomplete_at_3prime based on the matched region of the reference protein.')
    ap.add_argument('--genome', required=True, help='Sample genome FASTA (single-record FASTA preferred)')
    ap.add_argument('--gff', required=True, help='Input GFF3 after truncate-aware correction')
    ap.add_argument('--removed-log', required=True, help='TSV from truncate_aware_gff_v2.py listing gap-near-start annotations to remove')
    ap.add_argument('--reference-proteins', required=True, help='Reference proteins FASTA (e.g. 01_reference/all_refs_proteins.faa)')
    ap.add_argument('--out-gff', required=True, help='Output GFF3 with removed short artefacts excluded and rescued partial features added where supported')
    ap.add_argument('--recovered-log', required=True, help='TSV log of rescued partial annotations')
    ap.add_argument('--unrecovered-log', required=True, help='TSV log of removed genes that could not be recovered at the homology threshold')
    ap.add_argument('--min-query-coverage', type=float, default=0.50, help='Minimum fraction of the reference protein that must align to rescue a partial gene (default: 0.50)')
    ap.add_argument('--threads', type=int, default=1, help='Threads for tblastn (default: 1)')
    ap.add_argument('--evalue', type=float, default=1e-5, help='tblastn e-value threshold (default: 1e-5)')
    args = ap.parse_args()

    sample_name = find_sample_name(args.genome)
    removed_entries = parse_removed_log(args.removed_log, sample_name=sample_name)

    header, feats = parse_gff(args.gff)
    pruned_feats, actually_removed = remove_models_from_gff(feats, removed_entries)

    if not removed_entries:
        write_gff(header, pruned_feats, args.out_gff)
        with open(args.recovered_log, 'w', newline='') as out:
            writer = csv.DictWriter(out, fieldnames=['sample', 'gene', 'status'], delimiter='\t')
            writer.writeheader()
        with open(args.unrecovered_log, 'w', newline='') as out:
            writer = csv.DictWriter(out, fieldnames=['sample', 'gene', 'reason'], delimiter='\t')
            writer.writeheader()
        print('No removed annotations for this sample; copied GFF unchanged (after any direct removal checks).')
        return

    ref_by_gene = parse_reference_proteins(args.reference_proteins)
    existing_ids = set()
    for feat in pruned_feats:
        if 'ID' in feat['attrs']:
            existing_ids.add(feat['attrs']['ID'])

    recovered_rows = []
    unrecovered_rows = []
    recovered_features = []

    with tempfile.TemporaryDirectory(prefix='partial_recover_') as tmpdir:
        db_prefix = make_blast_db(args.genome, tmpdir)
        genes_to_rescue = []
        seen = set()
        for x in removed_entries:
            gene = x.get('gene')
            if gene and gene not in seen:
                genes_to_rescue.append(gene)
                seen.add(gene)

        for gene in genes_to_rescue:
            refs = ref_by_gene.get(gene, [])
            if not refs:
                unrecovered_rows.append({'sample': sample_name, 'gene': gene, 'reason': 'no_reference_protein_for_gene'})
                continue
            query_faa = os.path.join(tmpdir, f'{gene}.faa')
            out_tsv = os.path.join(tmpdir, f'{gene}.tblastn.tsv')
            write_query_fasta(refs, query_faa)
            try:
                run_tblastn(query_faa, db_prefix, out_tsv, threads=args.threads, evalue=args.evalue)
            except subprocess.CalledProcessError as e:
                unrecovered_rows.append({'sample': sample_name, 'gene': gene, 'reason': f'tblastn_failed:{e}'})
                continue
            hits = parse_tblastn_hits(out_tsv)
            best_chain = choose_best_partial_hit(hits, min_qcov=args.min_query_coverage)
            if best_chain is None:
                unrecovered_rows.append({'sample': sample_name, 'gene': gene, 'reason': f'no_chain_with_qcov>={args.min_query_coverage:.2f}'})
                continue
            ref_meta = next((rec for rec in refs if rec['header'] == best_chain['hits'][0]['qseqid']), refs[0])
            new_feats, log_row = build_partial_features(sample_name, gene, ref_meta, best_chain, existing_ids)
            recovered_features.extend(new_feats)
            recovered_rows.append({'sample': sample_name, **log_row, 'status': 'rescued_partial'})

    final_feats = pruned_feats + recovered_features
    write_gff(header, final_feats, args.out_gff)

    with open(args.recovered_log, 'w', newline='') as out:
        fieldnames = ['sample', 'gene', 'status', 'seqid', 'start', 'end', 'strand', 'blocks', 'incomplete_at_5prime', 'incomplete_at_3prime', 'qstart', 'qend', 'qlen', 'qcov', 'pident', 'ref_accession', 'ref_locus_tag', 'ref_product']
        writer = csv.DictWriter(out, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()
        for row in recovered_rows:
            writer.writerow(row)

    with open(args.unrecovered_log, 'w', newline='') as out:
        fieldnames = ['sample', 'gene', 'reason']
        writer = csv.DictWriter(out, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()
        for row in unrecovered_rows:
            writer.writerow(row)

    print('Done')
    print(f'Output GFF: {args.out_gff}')
    print(f'Recovered log: {args.recovered_log}')
    print(f'Unrecovered log: {args.unrecovered_log}')


if __name__ == '__main__':
    main()
