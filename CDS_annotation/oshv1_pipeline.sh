#!/usr/bin/env bash
set -euo pipefail

############################################################
# oshv1_pipeline.sh
#
# Modes:
#   download_refs  -> download OsHV-1 reference GBFF files from NCBI
#   annotate       -> annotate genomes from a single multifasta file
############################################################

usage() {
  cat <<'EOF'
Usage:
  ./oshv1_pipeline.sh <mode> [options]

Modes:
  download_refs   Download OsHV-1 reference GenBank files from NCBI into refs_gbff/
  annotate        Annotate genomes from a single multifasta file

------------------------------------------------------------
Mode: download_refs
------------------------------------------------------------
Required:
  -o   Output directory for downloaded reference GBFF files

Optional:
  -q   NCBI Entrez query
       Default:
       "Ostreid herpesvirus 1"[Organism] AND ("complete genome"[Title] OR "complete genome"[All Fields])

  -m   Maximum number of accessions to download (default: all)
  -f   Force redownload of existing files

Example:
  ./oshv1_pipeline.sh download_refs -o refs_gbff

------------------------------------------------------------
Mode: annotate
------------------------------------------------------------
Required:
  -i   Multifasta file containing genomes to annotate
  -r   Directory containing reference GenBank files (*.gb, *.gbk, *.gbff)
  -p   Primary reference accession for Liftoff (e.g. NC_005881.2)
  -o   Output directory

Optional:
  -t   Threads (default: 4)
  --min-ident   Minimum BLAST percent identity for adding missed CDS (default: 70)
  --min-qcov    Minimum BLAST query coverage (%) for adding missed CDS (default: 70)

Example:
  ./oshv1_pipeline.sh annotate \
    -i /group/sequencing/assembly/James/OsHV_annot/oshv1_refs.fasta \
    -r refs_gbff \
    -p NC_005881.2 \
    -o oshv1_annotation \
    -t 8
EOF
}

need_cmd() {
  command -v "$1" >/dev/null 2>&1 || {
    echo "ERROR: Missing dependency: $1" >&2
    exit 1
  }
}

write_split_multifasta_py() {
  local outfile="$1"
  cat > "$outfile" <<'PYEOF'
#!/usr/bin/env python3
import sys, os, re
from Bio import SeqIO

infile = sys.argv[1]
outdir = sys.argv[2]
os.makedirs(outdir, exist_ok=True)

seen = {}
manifest = os.path.join(outdir, "manifest.tsv")

def sanitize(s):
    s = s.strip()
    s = s.replace(" ", "_")
    s = re.sub(r'[^\w.\-]+', '_', s)
    return s or "unnamed_record"

with open(manifest, "w") as man:
    man.write("original_id\tsanitized_id\tlength\tdescription\n")
    n = 0
    for record in SeqIO.parse(infile, "fasta"):
        base = sanitize(record.id)
        if base in seen:
            seen[base] += 1
            sid = f"{base}_{seen[base]}"
        else:
            seen[base] = 1
            sid = base

        outfile = os.path.join(outdir, f"{sid}.fasta")
        with open(outfile, "w") as out:
            out.write(f">{sid}\n")
            seq = str(record.seq)
            for i in range(0, len(seq), 80):
                out.write(seq[i:i+80] + "\n")

        man.write(f"{record.id}\t{sid}\t{len(record.seq)}\t{record.description}\n")
        n += 1

if n == 0:
    sys.stderr.write("ERROR: No FASTA records found in input multifasta.\n")
    sys.exit(1)
PYEOF
  chmod +x "$outfile"
}

write_parse_refs_py() {
  local outfile="$1"
  cat > "$outfile" <<'PYEOF'
#!/usr/bin/env python3
import sys, os, re, glob
from Bio import SeqIO

refs_dir = sys.argv[1]
outdir = sys.argv[2]

os.makedirs(outdir, exist_ok=True)
ref_fna_dir = os.path.join(outdir, "ref_fastas")
ref_gff_dir = os.path.join(outdir, "ref_gffs")
os.makedirs(ref_fna_dir, exist_ok=True)
os.makedirs(ref_gff_dir, exist_ok=True)

all_prot_faa = os.path.join(outdir, "all_refs_proteins.faa")
manifest_tsv = os.path.join(outdir, "reference_manifest.tsv")

orf_regex = re.compile(r'\b(ORF[_\-\s]?[A-Za-z0-9.]+)\b', re.IGNORECASE)

def clean(s):
    if s is None:
        return "NA"
    s = str(s).strip()
    s = s.replace("\n", " ")
    s = re.sub(r'\s+', '_', s)
    s = re.sub(r'[^\w.\-|:+=]', '_', s)
    return s if s else "NA"

def extract_orf_name(feature):
    candidates = []
    for key in ["gene", "locus_tag", "product", "note", "label"]:
        if key in feature.qualifiers:
            candidates.extend(feature.qualifiers[key])
    joined = " ".join(candidates)
    m = orf_regex.search(joined)
    if m:
        val = m.group(1)
        val = re.sub(r'\s+', '', val.upper())
        val = val.replace("ORF_", "ORF").replace("ORF-", "ORF")
        return val
    gene = feature.qualifiers.get("gene", [""])[0]
    if gene.upper().startswith("ORF"):
        return re.sub(r'\s+', '', gene.upper())
    return None

def write_gff3_for_record(record, gff_path):
    with open(gff_path, "w") as gff:
        gff.write("##gff-version 3\n")
        for i, feat in enumerate(record.features, start=1):
            if feat.type != "CDS":
                continue
            if "translation" not in feat.qualifiers:
                continue

            seqid = record.id
            strand = "+" if feat.location.strand == 1 else "-"
            protein_id = feat.qualifiers.get("protein_id", [f"{record.id}_prot_{i}"])[0]
            locus_tag = feat.qualifiers.get("locus_tag", [f"{record.id}_cds_{i:04d}"])[0]
            gene = feat.qualifiers.get("gene", [locus_tag])[0]
            product = feat.qualifiers.get("product", ["hypothetical protein"])[0]
            orf_name = extract_orf_name(feat) or gene
            note = feat.qualifiers.get("note", [""])[0]

            if hasattr(feat.location, "parts"):
                parts = feat.location.parts
            else:
                parts = [feat.location]

            starts = [int(p.start) + 1 for p in parts]
            ends = [int(p.end) for p in parts]
            gene_start = min(starts)
            gene_end = max(ends)

            gene_id = f"{record.id}.gene.{i:04d}"
            mrna_id = f"{record.id}.mrna.{i:04d}"

            gene_attrs = {
                "ID": gene_id,
                "Name": clean(orf_name),
                "gene": clean(orf_name),
                "locus_tag": clean(locus_tag),
                "product": clean(product)
            }
            mrna_attrs = {
                "ID": mrna_id,
                "Parent": gene_id,
                "Name": clean(orf_name),
                "gene": clean(orf_name),
                "locus_tag": clean(locus_tag),
                "product": clean(product)
            }

            gff.write(
                f"{seqid}\tRefGBFF\tgene\t{gene_start}\t{gene_end}\t.\t{strand}\t.\t" +
                ";".join([f"{k}={v}" for k, v in gene_attrs.items()]) + "\n"
            )
            gff.write(
                f"{seqid}\tRefGBFF\tmRNA\t{gene_start}\t{gene_end}\t.\t{strand}\t.\t" +
                ";".join([f"{k}={v}" for k, v in mrna_attrs.items()]) + "\n"
            )

            for j, part in enumerate(parts, start=1):
                start = int(part.start) + 1
                end = int(part.end)
                phase = "0"
                cds_attrs = {
                    "ID": f"{mrna_id}.cds{j}",
                    "Parent": mrna_id,
                    "Name": clean(orf_name),
                    "gene": clean(orf_name),
                    "locus_tag": clean(locus_tag),
                    "product": clean(product),
                    "protein_id": clean(protein_id)
                }
                if note:
                    cds_attrs["note"] = clean(note)

                gff.write(
                    f"{seqid}\tRefGBFF\tCDS\t{start}\t{end}\t.\t{strand}\t{phase}\t" +
                    ";".join([f"{k}={v}" for k, v in cds_attrs.items()]) + "\n"
                )

gb_files = sorted(glob.glob(os.path.join(refs_dir, "*.gb")) +
                  glob.glob(os.path.join(refs_dir, "*.gbk")) +
                  glob.glob(os.path.join(refs_dir, "*.gbff")))

if not gb_files:
    sys.stderr.write("ERROR: No GenBank files found in reference directory.\n")
    sys.exit(1)

with open(all_prot_faa, "w") as prot_out, open(manifest_tsv, "w") as man:
    man.write("accession\tfeature_index\tlocus_tag\tgene\torf_name\tprotein_id\tproduct\n")

    for gb in gb_files:
        for record in SeqIO.parse(gb, "genbank"):
            accession = record.id
            fasta_out = os.path.join(ref_fna_dir, f"{accession}.fna")
            gff_out = os.path.join(ref_gff_dir, f"{accession}.gff3")

            with open(fasta_out, "w") as fna:
                fna.write(f">{record.id}\n")
                seq = str(record.seq)
                for i in range(0, len(seq), 80):
                    fna.write(seq[i:i+80] + "\n")

            write_gff3_for_record(record, gff_out)

            idx = 0
            for feat in record.features:
                if feat.type != "CDS":
                    continue
                if "translation" not in feat.qualifiers:
                    continue
                idx += 1

                locus_tag = feat.qualifiers.get("locus_tag", [f"{accession}_cds_{idx:04d}"])[0]
                gene = feat.qualifiers.get("gene", [locus_tag])[0]
                protein_id = feat.qualifiers.get("protein_id", [f"{accession}_prot_{idx:04d}"])[0]
                product = feat.qualifiers.get("product", ["hypothetical protein"])[0]
                orf_name = extract_orf_name(feat) or gene
                protein = feat.qualifiers["translation"][0].replace(" ", "").replace("\n", "")

                header = "|".join([
                    clean(accession),
                    clean(str(idx)),
                    clean(locus_tag),
                    clean(gene),
                    clean(orf_name),
                    clean(protein_id),
                    clean(product)
                ])
                prot_out.write(f">{header}\n")
                for i in range(0, len(protein), 80):
                    prot_out.write(protein[i:i+80] + "\n")

                man.write("\t".join([
                    accession, str(idx), locus_tag, gene, str(orf_name), protein_id, product
                ]) + "\n")
PYEOF
  chmod +x "$outfile"
}

write_merge_annotations_py() {
  local outfile="$1"
  cat > "$outfile" <<'PYEOF'
#!/usr/bin/env python3
import sys, re
from collections import defaultdict

liftoff_gff = sys.argv[1]
prodigal_gff = sys.argv[2]
blast_tsv = sys.argv[3]
output_gff = sys.argv[4]
summary_tsv = sys.argv[5]
min_ident = float(sys.argv[6])
min_qcov = float(sys.argv[7])

orf_regex = re.compile(r'\bORF[_\-\s]?[A-Za-z0-9.]+\b', re.IGNORECASE)

def parse_attrs(s):
    d = {}
    for item in s.strip().split(";"):
        if "=" in item:
            k, v = item.split("=", 1)
            d[k] = v
    return d

def attrs_to_str(d):
    ordered = []
    for k in ["ID", "Parent", "Name", "gene", "locus_tag", "product", "protein_id", "note"]:
        if k in d and d[k] not in [None, ""]:
            ordered.append(f"{k}={d[k]}")
    for k, v in d.items():
        if k not in {"ID", "Parent", "Name", "gene", "locus_tag", "product", "protein_id", "note"} and v not in [None, ""]:
            ordered.append(f"{k}={v}")
    return ";".join(ordered)

def extract_known_name(attrs):
    for key in ["gene", "Name", "locus_tag", "product", "note", "ID"]:
        if key in attrs:
            m = orf_regex.search(attrs[key])
            if m:
                x = m.group(0).upper()
                x = x.replace("ORF_", "ORF").replace("ORF-", "ORF").replace(" ", "")
                return x
    if "gene" in attrs:
        return attrs["gene"]
    if "Name" in attrs:
        return attrs["Name"]
    return None

def overlap_fraction(a_start, a_end, b_start, b_end):
    ov = max(0, min(a_end, b_end) - max(a_start, b_start) + 1)
    if ov == 0:
        return 0.0
    a_len = a_end - a_start + 1
    return ov / a_len

header_lines = []
features = []
existing_names = set()
existing_intervals = defaultdict(list)

with open(liftoff_gff) as fh:
    for line in fh:
        if line.startswith("#"):
            header_lines.append(line.rstrip("\n"))
            continue
        parts = line.rstrip("\n").split("\t")
        if len(parts) != 9:
            continue
        seqid, source, ftype, start, end, score, strand, phase, attr = parts
        start, end = int(start), int(end)
        attrs = parse_attrs(attr)
        name = extract_known_name(attrs)
        if name:
            existing_names.add(name)
        if ftype in {"gene", "mRNA", "CDS"}:
            existing_intervals[seqid].append((start, end))
        features.append((seqid, source, ftype, start, end, score, strand, phase, attrs))

prodigal = {}
with open(prodigal_gff) as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        parts = line.rstrip("\n").split("\t")
        if len(parts) != 9:
            continue
        seqid, source, ftype, start, end, score, strand, phase, attr = parts
        if ftype != "CDS":
            continue
        attrs = parse_attrs(attr)
        prod_id = attrs.get("ID")
        if not prod_id:
            continue
        prodigal[prod_id] = {
            "seqid": seqid,
            "start": int(start),
            "end": int(end),
            "strand": strand,
            "phase": phase if phase != "." else "0",
            "attrs": attrs
        }

best = {}
with open(blast_tsv) as fh:
    for line in fh:
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 15:
            continue
        qseqid, sseqid, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore, qlen, slen, qcovs = parts[:15]
        pident = float(pident)
        bitscore = float(bitscore)
        qcovs = float(qcovs)

        prev = best.get(qseqid)
        if prev is None or bitscore > prev["bitscore"]:
            fields = sseqid.split("|")
            ref = {
                "accession": fields[0] if len(fields) > 0 else "NA",
                "feature_index": fields[1] if len(fields) > 1 else "NA",
                "locus_tag": fields[2] if len(fields) > 2 else "NA",
                "gene": fields[3] if len(fields) > 3 else "NA",
                "orf_name": fields[4] if len(fields) > 4 else "NA",
                "protein_id": fields[5] if len(fields) > 5 else "NA",
                "product": fields[6] if len(fields) > 6 else "hypothetical_protein",
            }
            best[qseqid] = {
                "pident": pident,
                "bitscore": bitscore,
                "qcovs": qcovs,
                "ref": ref,
                "evalue": evalue
            }

added = []
gene_counter = 1

for qid, info in best.items():
    if qid not in prodigal:
        continue
    if info["pident"] < min_ident or info["qcovs"] < min_qcov:
        continue

    hit = info["ref"]
    proposed_name = hit["orf_name"] if hit["orf_name"] not in {"NA", ""} else hit["gene"]
    proposed_name = proposed_name.replace("ORF_", "ORF").replace("ORF-", "ORF").replace(" ", "")

    if proposed_name in existing_names:
        continue

    p = prodigal[qid]
    seqid, start, end = p["seqid"], p["start"], p["end"]

    overlaps = [overlap_fraction(start, end, s, e) for (s, e) in existing_intervals.get(seqid, [])]
    if overlaps and max(overlaps) >= 0.50:
        continue

    gene_id = f"added_homology_gene_{gene_counter:04d}"
    mrna_id = f"added_homology_mrna_{gene_counter:04d}"
    note = (
        f"inferred_from_homology;best_hit={hit['accession']}|{hit['locus_tag']}|{hit['orf_name']};"
        f"pident={info['pident']:.2f};qcov={info['qcovs']:.2f};evalue={info['evalue']}"
    )

    gene_attrs = {
        "ID": gene_id,
        "Name": proposed_name,
        "gene": proposed_name,
        "locus_tag": qid,
        "product": hit["product"]
    }
    mrna_attrs = {
        "ID": mrna_id,
        "Parent": gene_id,
        "Name": proposed_name,
        "gene": proposed_name,
        "locus_tag": qid,
        "product": hit["product"]
    }
    cds_attrs = {
        "ID": f"{mrna_id}.cds1",
        "Parent": mrna_id,
        "Name": proposed_name,
        "gene": proposed_name,
        "locus_tag": qid,
        "product": hit["product"],
        "protein_id": qid,
        "note": note
    }

    features.append((seqid, "Homology+Prodigal", "gene", start, end, ".", p["strand"], ".", gene_attrs))
    features.append((seqid, "Homology+Prodigal", "mRNA", start, end, ".", p["strand"], ".", mrna_attrs))
    features.append((seqid, "Homology+Prodigal", "CDS", start, end, ".", p["strand"], "0", cds_attrs))

    existing_names.add(proposed_name)
    existing_intervals[seqid].append((start, end))

    added.append({
        "seqid": seqid,
        "start": start,
        "end": end,
        "strand": p["strand"],
        "assigned_name": proposed_name,
        "ref_accession": hit["accession"],
        "ref_locus_tag": hit["locus_tag"],
        "ref_orf_name": hit["orf_name"],
        "product": hit["product"],
        "pident": info["pident"],
        "qcovs": info["qcovs"],
        "evalue": info["evalue"]
    })
    gene_counter += 1

type_order = {"gene": 0, "mRNA": 1, "CDS": 2}
features_sorted = sorted(features, key=lambda x: (x[0], x[3], x[4], type_order.get(x[2], 9)))

with open(output_gff, "w") as out:
    if header_lines:
        for h in header_lines:
            out.write(h + "\n")
    else:
        out.write("##gff-version 3\n")
    for seqid, source, ftype, start, end, score, strand, phase, attrs in features_sorted:
        out.write("\t".join([
            seqid, source, ftype, str(start), str(end), score, strand, phase, attrs_to_str(attrs)
        ]) + "\n")

with open(summary_tsv, "w") as out:
    out.write("seqid\tstart\tend\tstrand\tassigned_name\tref_accession\tref_locus_tag\tref_orf_name\tproduct\tpident\tqcovs\tevalue\n")
    for row in added:
        out.write("\t".join([
            str(row["seqid"]), str(row["start"]), str(row["end"]), str(row["strand"]),
            str(row["assigned_name"]), str(row["ref_accession"]), str(row["ref_locus_tag"]),
            str(row["ref_orf_name"]), str(row["product"]),
            f"{row['pident']:.2f}", f"{row['qcovs']:.2f}", str(row["evalue"])
        ]) + "\n")
PYEOF
  chmod +x "$outfile"
}

download_refs_mode() {
  local OUTDIR=""
  local QUERY='"Ostreid herpesvirus 1"[Organism] AND ("complete genome"[Title] OR "complete genome"[All Fields])'
  local MAXN=""
  local FORCE=0

  while [[ $# -gt 0 ]]; do
    case "$1" in
      -o) OUTDIR="$2"; shift 2 ;;
      -q) QUERY="$2"; shift 2 ;;
      -m) MAXN="$2"; shift 2 ;;
      -f) FORCE=1; shift ;;
      -h|--help) usage; exit 0 ;;
      *) echo "ERROR: Unknown option for download_refs: $1" >&2; usage; exit 1 ;;
    esac
  done

  if [[ -z "$OUTDIR" ]]; then
    echo "ERROR: download_refs requires -o" >&2
    usage
    exit 1
  fi

  need_cmd esearch
  need_cmd efetch

  mkdir -p "$OUTDIR"
  mkdir -p "$OUTDIR/_logs"

  local acc_file="$OUTDIR/_logs/ncbi_accessions.txt"

  echo "[INFO] Querying NCBI..."
  echo "[INFO] Query: $QUERY"

  esearch -db nucleotide -query "$QUERY" | efetch -format acc > "$acc_file"

  if [[ ! -s "$acc_file" ]]; then
    echo "ERROR: No accessions returned by NCBI query." >&2
    exit 1
  fi

  if [[ -n "$MAXN" ]]; then
    head -n "$MAXN" "$acc_file" > "${acc_file}.tmp"
    mv "${acc_file}.tmp" "$acc_file"
  fi

  echo "[INFO] Accessions:"
  cat "$acc_file"
  echo

  local n=0
  while read -r acc; do
    [[ -z "$acc" ]] && continue
    local outfile="$OUTDIR/${acc}.gbff"

    if [[ -s "$outfile" && "$FORCE" -eq 0 ]]; then
      echo "[SKIP] $outfile exists"
      continue
    fi

    echo "[INFO] Downloading $acc"
    efetch -db nucleotide -id "$acc" -format gbwithparts > "$outfile"

    if [[ ! -s "$outfile" ]]; then
      echo "ERROR: Failed to download $acc" >&2
      exit 1
    fi

    n=$((n+1))
  done < "$acc_file"

  echo
  echo "[DONE] References downloaded into: $OUTDIR"
  echo "[INFO] Newly downloaded this run: $n"
}

annotate_mode() {
  local MULTIFASTA=""
  local REFS_DIR=""
  local PRIMARY_REF=""
  local OUTDIR=""
  local THREADS=4
  local MIN_IDENT=70
  local MIN_QCOV=70

  while [[ $# -gt 0 ]]; do
    case "$1" in
      -i) MULTIFASTA="$2"; shift 2 ;;
      -r) REFS_DIR="$2"; shift 2 ;;
      -p) PRIMARY_REF="$2"; shift 2 ;;
      -o) OUTDIR="$2"; shift 2 ;;
      -t) THREADS="$2"; shift 2 ;;
      --min-ident) MIN_IDENT="$2"; shift 2 ;;
      --min-qcov) MIN_QCOV="$2"; shift 2 ;;
      -h|--help) usage; exit 0 ;;
      *) echo "ERROR: Unknown option for annotate: $1" >&2; usage; exit 1 ;;
    esac
  done

  if [[ -z "$MULTIFASTA" || -z "$REFS_DIR" || -z "$PRIMARY_REF" || -z "$OUTDIR" ]]; then
    echo "ERROR: annotate requires -i, -r, -p, -o" >&2
    usage
    exit 1
  fi

  [[ -s "$MULTIFASTA" ]] || { echo "ERROR: Missing input multifasta: $MULTIFASTA" >&2; exit 1; }
  [[ -d "$REFS_DIR" ]] || { echo "ERROR: Missing refs dir: $REFS_DIR" >&2; exit 1; }

  need_cmd python
  need_cmd prodigal
  need_cmd makeblastdb
  need_cmd blastp
  need_cmd liftoff
  need_cmd minimap2

  mkdir -p "$OUTDIR"/{00_logs,00_split_genomes,01_reference,02_predictions,03_liftoff,04_blast,05_final,06_qc,tmp}

  write_split_multifasta_py "$OUTDIR/tmp/split_multifasta.py"
  write_parse_refs_py "$OUTDIR/tmp/parse_refs.py"
  write_merge_annotations_py "$OUTDIR/tmp/merge_annotations.py"

  echo "[INFO] Splitting multifasta..."
  python "$OUTDIR/tmp/split_multifasta.py" "$MULTIFASTA" "$OUTDIR/00_split_genomes"

  echo "[INFO] Parsing reference GBFF files..."
  python "$OUTDIR/tmp/parse_refs.py" "$REFS_DIR" "$OUTDIR/01_reference"

  [[ -s "$OUTDIR/01_reference/all_refs_proteins.faa" ]] || {
    echo "ERROR: Failed to build reference protein FASTA" >&2
    exit 1
  }

  local PRIMARY_REF_FASTA="$OUTDIR/01_reference/ref_fastas/${PRIMARY_REF}.fna"
  local PRIMARY_REF_GFF="$OUTDIR/01_reference/ref_gffs/${PRIMARY_REF}.gff3"

  [[ -s "$PRIMARY_REF_FASTA" ]] || {
    echo "ERROR: Primary reference FASTA not found: $PRIMARY_REF_FASTA" >&2
    echo "Check the exact accession used with -p" >&2
    exit 1
  }

  [[ -s "$PRIMARY_REF_GFF" ]] || {
    echo "ERROR: Primary reference GFF not found: $PRIMARY_REF_GFF" >&2
    exit 1
  }

  echo "[INFO] Building BLAST database..."
  makeblastdb \
    -in "$OUTDIR/01_reference/all_refs_proteins.faa" \
    -dbtype prot \
    -out "$OUTDIR/01_reference/oshv1_ref_proteins_db" \
    > "$OUTDIR/00_logs/makeblastdb.log" 2>&1

  shopt -s nullglob
  local GENOME_FILES=("$OUTDIR/00_split_genomes"/*.fasta)
  shopt -u nullglob

  [[ ${#GENOME_FILES[@]} -gt 0 ]] || {
    echo "ERROR: No split genome FASTA files found" >&2
    exit 1
  }

  echo "[INFO] Found ${#GENOME_FILES[@]} genome(s) to annotate"

  local genome
  for genome in "${GENOME_FILES[@]}"; do
    local sample
    sample=$(basename "$genome")
    sample="${sample%.fasta}"

    echo "[INFO] Processing $sample"

    mkdir -p "$OUTDIR/02_predictions/$sample" \
             "$OUTDIR/03_liftoff/$sample" \
             "$OUTDIR/04_blast/$sample" \
             "$OUTDIR/05_final/$sample" \
             "$OUTDIR/06_qc/$sample"

    local prodigal_gff="$OUTDIR/02_predictions/$sample/${sample}.prodigal.gff"
    local prodigal_faa="$OUTDIR/02_predictions/$sample/${sample}.prodigal.faa"
    local prodigal_fna="$OUTDIR/02_predictions/$sample/${sample}.prodigal.fna"

    prodigal \
      -i "$genome" \
      -o "$prodigal_gff" \
      -a "$prodigal_faa" \
      -d "$prodigal_fna" \
      -f gff \
      -p meta \
      > "$OUTDIR/00_logs/${sample}.prodigal.log" 2>&1

    local blast_tsv="$OUTDIR/04_blast/$sample/${sample}.vs_refs.blastp.tsv"

    blastp \
      -query "$prodigal_faa" \
      -db "$OUTDIR/01_reference/oshv1_ref_proteins_db" \
      -out "$blast_tsv" \
      -evalue 1e-5 \
      -max_target_seqs 5 \
      -num_threads "$THREADS" \
      -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovs" \
      > "$OUTDIR/00_logs/${sample}.blastp.log" 2>&1

    local liftoff_gff="$OUTDIR/03_liftoff/$sample/${sample}.liftoff.gff3"
    local liftoff_unmapped="$OUTDIR/03_liftoff/$sample/${sample}.liftoff.unmapped.txt"

    liftoff "$genome" "$PRIMARY_REF_FASTA" \
      -g "$PRIMARY_REF_GFF" \
      -o "$liftoff_gff" \
      -u "$liftoff_unmapped" \
      -p "$THREADS" \
      -copies \
      > "$OUTDIR/00_logs/${sample}.liftoff.log" 2>&1

    local final_gff="$OUTDIR/05_final/$sample/${sample}.final.gff3"
    local added_summary="$OUTDIR/06_qc/$sample/${sample}.homology_added.tsv"

    python "$OUTDIR/tmp/merge_annotations.py" \
      "$liftoff_gff" \
      "$prodigal_gff" \
      "$blast_tsv" \
      "$final_gff" \
      "$added_summary" \
      "$MIN_IDENT" \
      "$MIN_QCOV"

    grep -P '\tCDS\t' "$liftoff_gff" | wc -l > "$OUTDIR/06_qc/$sample/${sample}.liftoff_cds_count.txt" || true
    grep -P '\tCDS\t' "$final_gff" | wc -l > "$OUTDIR/06_qc/$sample/${sample}.final_cds_count.txt" || true

    awk -F'\t' 'BEGIN{OFS="\t"} !/^#/ && $3=="CDS" {print $1,$4,$5,$7,$9}' "$final_gff" \
      > "$OUTDIR/06_qc/$sample/${sample}.final_cds_coordinates.tsv"
  done

  echo
  echo "[DONE] Annotation complete"
  echo "  Split genomes:     $OUTDIR/00_split_genomes/"
  echo "  Reference parse:   $OUTDIR/01_reference/"
  echo "  Final GFF3 files:  $OUTDIR/05_final/*/*.final.gff3"
  echo "  QC summaries:      $OUTDIR/06_qc/"
}

if [[ $# -lt 1 ]]; then
  usage
  exit 1
fi

MODE="$1"
shift

case "$MODE" in
  download_refs) download_refs_mode "$@" ;;
  annotate) annotate_mode "$@" ;;
  -h|--help) usage ;;
  *)
    echo "ERROR: Unknown mode: $MODE" >&2
    usage
    exit 1
    ;;
esac
