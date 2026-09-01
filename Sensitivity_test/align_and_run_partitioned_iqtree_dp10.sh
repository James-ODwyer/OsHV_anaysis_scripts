#!/bin/bash
#SBATCH --nodes=1
#SBATCH --job-name="MAFFT_iqtree"
#SBATCH --account="aggstrategic1"
#SBATCH --partition="batch"
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=48GB
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=james.o\'dwyer@agriculture.vic.gov.au
#SBATCH --time=72:0:00


set -euo pipefail

module load MAFFT/7.526-GCC-13.3.0-with-extensions

module load IQ-TREE/2.3.6-gompi-2024a

module load Python/3.12.3-GCCcore-13.3.0

module load Python-bundle-PyPI


eval "$(conda shell.bash hook)"
conda activate /home/vidh72t/.conda/envs/oshv1_annot


# Align each retained OsHV-1 CDS independently with MAFFT, concatenate the
# alignments, construct an IQ-TREE NEXUS partition file, and infer a
# partitioned maximum-likelihood tree with an independently selected
# substitution model for each CDS partition.
#
# The MAFFT parameters reproduce those described in the manuscript:
#   --op 2 --ep 0.12 --maxiterate 1000 --retree 150
#
# IQ-TREE uses -p, which applies a partitioned analysis with a separate
# substitution model selected for each gene while retaining proportional
# branch lengths across partitions. Set BRANCH_MODE="-sp" below if fully
# independent branch lengths are specifically required instead.



#INPUT_DIR="/group/sequencing/assembly/James/OsHV_annot/Align_raws/dp10_complete_CDS_sets/per_gene_at_least_10_samples"
#INPUT_DIR="/group/sequencing/assembly/James/OsHV_annot/Align_raws/dp10_complete_CDS_sets/per_gene_at_least_10_samples_added_SRA"
INPUT_DIR="/group/sequencing/assembly/James/OsHV_annot/Align_raws/dp10_complete_CDS_sets/per_gene_at_least_10_samples_added_SRA_nodups"


OUTPUT_DIR="/group/sequencing/assembly/James/OsHV_annot/Align_raws/dp10_complete_CDS_sets/mafft_iqtree_partitioned_dp10_added_SRA_nodups"

THREADS="${THREADS:-AUTO}"
BOOTSTRAPS="${BOOTSTRAPS:-2000}"
SH_ALRT="${SH_ALRT:-2000}"
SEED="${SEED:-12345}"

# -p: separate best-fit model per partition, proportional branch lengths.
# -sp: separate best-fit model and fully independent branch lengths.
BRANCH_MODE="${BRANCH_MODE:--p}"

# IQ-TREE executable may be called iqtree2 or iqtree depending on installation.
if command -v iqtree2 >/dev/null 2>&1; then
    IQTREE_BIN="iqtree2"
elif command -v iqtree >/dev/null 2>&1; then
    IQTREE_BIN="iqtree"
else
    echo "ERROR: neither iqtree2 nor iqtree is available in PATH." >&2
    exit 1
fi

if ! command -v mafft >/dev/null 2>&1; then
    echo "ERROR: mafft is not available in PATH." >&2
    exit 1
fi

if [[ ! -d "$INPUT_DIR" ]]; then
    echo "ERROR: input directory does not exist: $INPUT_DIR" >&2
    exit 1
fi

case "$BRANCH_MODE" in
    -p|-sp) ;;
    *)
        echo "ERROR: BRANCH_MODE must be -p or -sp, not: $BRANCH_MODE" >&2
        exit 1
        ;;
esac

ALIGN_DIR="$OUTPUT_DIR/individual_gene_alignments"
LOG_DIR="$OUTPUT_DIR/logs"
CONCAT_DIR="$OUTPUT_DIR/concatenated"
mkdir -p "$ALIGN_DIR" "$LOG_DIR" "$CONCAT_DIR"
# Prevent stale alignments from a previous run from entering the concatenation.
rm -f "$ALIGN_DIR"/*.aligned.fasta "$LOG_DIR"/*.mafft.log

MANIFEST="$OUTPUT_DIR/gene_alignment_manifest.tsv"
PARTITIONS="$CONCAT_DIR/partitions.nex"
CONCAT_FASTA="$CONCAT_DIR/concatenated_CDS_alignment.fasta"
PREFIX="$CONCAT_DIR/OsHV_dp10_partitioned"

printf 'gene\tinput_fasta\talignment_fasta\tn_sequences\talignment_length\tstart\tend\n' > "$MANIFEST"

mapfile -d '' FASTAS < <(
    find "$INPUT_DIR" -maxdepth 1 -type f \
        \( -iname '*.fa' -o -iname '*.fas' -o -iname '*.fasta' -o -iname '*.fna' \) \
        -print0 | sort -zV
)

if (( ${#FASTAS[@]} == 0 )); then
    echo "ERROR: no FASTA files were found in $INPUT_DIR" >&2
    exit 1
fi

echo "Found ${#FASTAS[@]} gene FASTA files."

# Align every gene independently using the manuscript parameters. CDS
# sequences are already extracted in coding orientation, so sequence IDs and
# orientations are retained unchanged.
for fasta in "${FASTAS[@]}"; do
    filename="$(basename "$fasta")"
    gene="${filename%.*}"
    safe_gene="$(printf '%s' "$gene" | sed 's/[^A-Za-z0-9_.-]/_/g')"
    aln="$ALIGN_DIR/${safe_gene}.aligned.fasta"
    log="$LOG_DIR/${safe_gene}.mafft.log"

    echo "Aligning $gene"
    mafft \
        --op 2 \
        --ep 0.12 \
        --maxiterate 1000 \
        --retree 150 \
        --thread -1 \
        "$fasta" > "$aln" 2> "$log"

done

# Concatenate aligned loci by sequence ID. Missing loci are represented by N characters (unknown nucleotides), not alignment gaps.
# This is required for the >=10-sample dataset, where some genes are absent from
# some study samples. Biopython is used only for deterministic FASTA parsing and
# writing; no sequence transformation is performed.
python3 - "$ALIGN_DIR" "$CONCAT_FASTA" "$PARTITIONS" "$MANIFEST" "$INPUT_DIR" <<'PY'
from pathlib import Path
from collections import OrderedDict
import csv
import re
import sys

try:
    from Bio import SeqIO
except ImportError:
    raise SystemExit("ERROR: Biopython is required: python3 -m pip install biopython")

align_dir = Path(sys.argv[1])
concat_fasta = Path(sys.argv[2])
partition_file = Path(sys.argv[3])
manifest_file = Path(sys.argv[4])
input_dir = Path(sys.argv[5])

alignment_files = sorted(
    align_dir.glob("*.aligned.fasta"),
    key=lambda p: [int(x) if x.isdigit() else x.lower() for x in re.split(r"(\d+)", p.name)],
)
if not alignment_files:
    raise SystemExit("ERROR: no MAFFT alignments were generated")

loci = []
all_ids = OrderedDict()
for aln_path in alignment_files:
    records = list(SeqIO.parse(str(aln_path), "fasta"))
    if not records:
        raise SystemExit(f"ERROR: empty alignment: {aln_path}")

    lengths = {len(record.seq) for record in records}
    if len(lengths) != 1:
        raise SystemExit(f"ERROR: sequences have unequal lengths after alignment: {aln_path}")

    # Use the sample/reference identifier preceding the first pipe. The prior
    # extraction script writes headers as sample|gene|metadata.
    seqs = OrderedDict()
    for record in records:
        sequence_id = record.id.split("|", 1)[0]
        if sequence_id in seqs:
            raise SystemExit(
                f"ERROR: duplicate normalized sequence ID '{sequence_id}' in {aln_path}. "
                "Inspect repeated ORF copies or FASTA headers."
            )
        seqs[sequence_id] = str(record.seq).upper()
        all_ids.setdefault(sequence_id, None)

    gene = aln_path.name.removesuffix(".aligned.fasta")
    loci.append((gene, aln_path, seqs, lengths.pop()))

concatenated = OrderedDict((seq_id, []) for seq_id in all_ids)
partition_rows = []
position = 1

for gene, aln_path, seqs, aln_length in loci:
    start = position
    end = position + aln_length - 1
    for seq_id in concatenated:
        concatenated[seq_id].append(seqs.get(seq_id, "N" * aln_length))
    partition_rows.append((gene, aln_path, len(seqs), aln_length, start, end))
    position = end + 1

with concat_fasta.open("w") as handle:
    for seq_id, pieces in concatenated.items():
        sequence = "".join(pieces)
        handle.write(f">{seq_id}\n")
        for i in range(0, len(sequence), 80):
            handle.write(sequence[i:i+80] + "\n")

# NEXUS SETS format accepted by IQ-TREE. DNA declares nucleotide partitions.
with partition_file.open("w") as handle:
    handle.write("#nexus\n")
    handle.write("begin sets;\n")
    for gene, _, _, _, start, end in partition_rows:
        partition_name = re.sub(r"[^A-Za-z0-9_.-]", "_", gene)
        handle.write(f"    charset {partition_name} = {start}-{end};\n")
    handle.write("end;\n")

# Add coordinates and alignment statistics to the manifest already created by bash.
with manifest_file.open("a", newline="") as handle:
    writer = csv.writer(handle, delimiter="\t")
    for gene, aln_path, count, length, start, end in partition_rows:
        candidates = list(input_dir.glob(gene + ".*"))
        input_path = str(candidates[0]) if candidates else ""
        writer.writerow([gene, input_path, aln_path, count, length, start, end])

expected_length = position - 1
observed_lengths = {sum(len(piece) for piece in pieces) for pieces in concatenated.values()}
if observed_lengths != {expected_length}:
    raise SystemExit(
        f"ERROR: concatenation length mismatch: expected {expected_length}, observed {observed_lengths}"
    )

print(f"Concatenated {len(loci)} loci for {len(concatenated)} sequence IDs.")
print(f"Final alignment length: {expected_length} nt.")
PY

# IQ-TREE partitioned analysis. MFP+MERGE permits selection of a separate
# best-fitting model for each partition and optionally merges statistically
# similar partitions. Replace with -m MFP if every input gene must remain a
# strictly separate model partition with no partition merging.
#
# Here we use -m MFP to preserve the user's selected genes as separate
# partitions and select the best model independently for each one.
echo "Running partitioned IQ-TREE analysis with $BRANCH_MODE"
"$IQTREE_BIN" \
    -s "$CONCAT_FASTA" \
    "$BRANCH_MODE" "$PARTITIONS" \
    -m MFP \
    -B "$BOOTSTRAPS" \
    --alrt "$SH_ALRT" \
    -T "$THREADS" \
    -seed "$SEED" \
    --prefix "$PREFIX" \
    --redo \
    2>&1 | tee "$LOG_DIR/iqtree_partitioned.log"

echo
echo "Analysis complete."
echo "Individual alignments: $ALIGN_DIR"
echo "Concatenated alignment: $CONCAT_FASTA"
echo "Partition file: $PARTITIONS"
echo "Input/alignment manifest: $MANIFEST"
echo "IQ-TREE prefix: $PREFIX"
echo "Main treefile: ${PREFIX}.treefile"
echo "Selected per-partition models: ${PREFIX}.best_scheme.nex (when emitted by IQ-TREE)"
 
