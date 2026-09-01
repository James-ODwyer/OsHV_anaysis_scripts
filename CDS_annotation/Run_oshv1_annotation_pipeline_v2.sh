#!/bin/bash
#SBATCH --nodes=1
#SBATCH --job-name="annotate_genomes"
#SBATCH --account="acct"
#SBATCH --partition="batch"
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=90GB
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=user@usermail.com
#SBATCH --time=6:0:00

set -euo pipefail

eval "$(conda shell.bash hook)"
conda activate /*/.conda/envs/oshv1_annot

############################################################
# GLOBAL VARIABLES
############################################################

BASE_DIR="/*/James/OsHV_annot"
INPUT_FASTA="${BASE_DIR}/OsHV_all_seq_genomes.fasta"
REFS_DIR="${BASE_DIR}/refs_gbff"

THREADS=24
MIN_IDENT=75
MIN_QCOV=75

PIPELINE="${BASE_DIR}/oshv1_pipeline.sh"
TRUNCATE_AWARE="${BASE_DIR}/truncate_aware_gff_v2.py"
COMPARE_MULTI="${BASE_DIR}/compare_oshv1_multi_runs.py"
MERGE_SAFE="${BASE_DIR}/merge_multiref_gffs_safe.py"
RECOVER_PARTIAL="${BASE_DIR}/recover_gap_partial_genes_v2.py"
PREMSTOP_SUMMARY="${BASE_DIR}/summarize_premature_stops_v2.py"

COMPARE_OUT="${BASE_DIR}/oshv1_three_run_comparison"
MERGED_OUT="${BASE_DIR}/oshv1_annotation_merged_MF509813plus"
PREMSTOP_OUT="${BASE_DIR}/oshv1_premature_stop_summary"

#
# RUN 1 VARIABLES: MF509813 PRIMARY

RUN1_LABEL="MF509813"
RUN1_ACC="MF509813.1"
RUN1_OUT="${BASE_DIR}/oshv1_annotation_MF509813_primary"
RUN1_REF_GFF="${RUN1_OUT}/01_reference/ref_gffs/${RUN1_ACC}.gff3"
RUN1_REF_PROTEINS="${RUN1_OUT}/01_reference/all_refs_proteins.faa"

###### RUN 2 VARIABLES: NC_005881 PRIMARY


RUN2_LABEL="NC_005881"
RUN2_ACC="NC_005881.2"
RUN2_OUT="${BASE_DIR}/oshv1_annotation_NC_005881"
RUN2_REF_GFF="${RUN2_OUT}/01_reference/ref_gffs/${RUN2_ACC}.gff3"
RUN2_REF_PROTEINS="${RUN2_OUT}/01_reference/all_refs_proteins.faa"


# RUN 3 VARIABLES: KY242785 PRIMARY

RUN3_LABEL="KY242785"
RUN3_ACC="KY242785.1"
RUN3_OUT="${BASE_DIR}/oshv1_annotation_KY242785"
RUN3_REF_GFF="${RUN3_OUT}/01_reference/ref_gffs/${RUN3_ACC}.gff3"
RUN3_REF_PROTEINS="${RUN3_OUT}/01_reference/all_refs_proteins.faa"

finalize_run() {
  local RUN_LABEL="$1"
  local RUN_ACC="$2"
  local RUN_OUT="$3"
  local RUN_REF_GFF="$4"
  local RUN_REF_PROTEINS="$5"

  echo "[INFO] =========================================================="
  echo "[INFO] Running annotation for ${RUN_LABEL} (${RUN_ACC})"
  echo "[INFO] Output directory: ${RUN_OUT}"
  echo "[INFO] =========================================================="

  ############################################################
  # STEP 1: MAIN ANNOTATION
  ############################################################
  "${PIPELINE}" annotate \
    -i "${INPUT_FASTA}" \
    -r "${REFS_DIR}" \
    -p "${RUN_ACC}" \
    -o "${RUN_OUT}" \
    --min-ident "${MIN_IDENT}" \
    --min-qcov "${MIN_QCOV}" \
    -t "${THREADS}"

  [[ -s "${RUN_REF_GFF}" ]] || { echo "[ERROR] Missing reference GFF for ${RUN_LABEL}: ${RUN_REF_GFF}" >&2; exit 1; }
  [[ -s "${RUN_REF_PROTEINS}" ]] || { echo "[ERROR] Missing reference proteins for ${RUN_LABEL}: ${RUN_REF_PROTEINS}" >&2; exit 1; }

  # STEP 2: TRUNCATION-AWARE CORRECTION + GAP-CAUSED REMOVAL

  echo "[INFO] Applying truncate-aware correction and start-gap removal for ${RUN_LABEL}"

  for gff in "${RUN_OUT}"/05_final/*/*.final.gff3; do
    sample=$(basename "${gff}" .final.gff3)
    genome="${RUN_OUT}/00_split_genomes/${sample}.fasta"
    outgff="$(dirname "${gff}")/${sample}.final.truncate_aware.gff3"
    summary="${RUN_OUT}/06_qc/${sample}/${sample}.truncate_aware.tsv"
    removed_log="${RUN_OUT}/06_qc/${sample}/${sample}.removed_start_gap_annotations.tsv"

    [[ -s "${genome}" ]] || { echo "[WARN] Missing genome FASTA for ${sample}, skipping truncate-aware step"; continue; }
    mkdir -p "$(dirname "${summary}")"

    python "${TRUNCATE_AWARE}" \
      --genome "${genome}" \
      --gff "${gff}" \
      --out-gff "${outgff}" \
      --summary "${summary}" \
      --reference-gff "${RUN_REF_GFF}" \
      --removed-log "${removed_log}" \
      --gap-within-bp 3 \
      --min-fraction-of-reference 0.90

    mv "${outgff}" "${gff}"
  done

  echo "[INFO] Recovering partial gap genes for ${RUN_LABEL}"

  for gff in "${RUN_OUT}"/05_final/*/*.final.gff3; do
    sample=$(basename "${gff}" .final.gff3)

    genome="${RUN_OUT}/00_split_genomes/${sample}.fasta"
    removed_log="${RUN_OUT}/06_qc/${sample}/${sample}.removed_start_gap_annotations.tsv"

    recovered_gff="$(dirname "${gff}")/${sample}.final.partial_rescued.gff3"
    recovered_log="${RUN_OUT}/06_qc/${sample}/${sample}.recovered_gap_partials.tsv"
    unrecovered_log="${RUN_OUT}/06_qc/${sample}/${sample}.unrecovered_gap_partials.tsv"

    [[ -s "${genome}" ]] || { echo "[WARN] Missing genome FASTA for ${sample}, skipping recovery"; continue; }
    [[ -s "${removed_log}" ]] || { echo "[INFO] No removed-start-gap log for ${sample}, skipping recovery"; continue; }

    python "${RECOVER_PARTIAL}" \
      --genome "${genome}" \
      --gff "${gff}" \
      --removed-log "${removed_log}" \
      --reference-proteins "${RUN_REF_PROTEINS}" \
      --out-gff "${recovered_gff}" \
      --recovered-log "${recovered_log}" \
      --unrecovered-log "${unrecovered_log}" \
      --min-query-coverage 0.50 \
      --threads "${THREADS}"

    mv "${recovered_gff}" "${gff}"
  done

  echo "[INFO] Finalised ${RUN_LABEL}"
}


# Run all three references 


finalize_run "${RUN1_LABEL}" "${RUN1_ACC}" "${RUN1_OUT}" "${RUN1_REF_GFF}" "${RUN1_REF_PROTEINS}"
finalize_run "${RUN2_LABEL}" "${RUN2_ACC}" "${RUN2_OUT}" "${RUN2_REF_GFF}" "${RUN2_REF_PROTEINS}"
finalize_run "${RUN3_LABEL}" "${RUN3_ACC}" "${RUN3_OUT}" "${RUN3_REF_GFF}" "${RUN3_REF_PROTEINS}"


echo "[INFO] =========================================================="
echo "[INFO] Comparing finalised runs"
echo "[INFO] =========================================================="

python "${COMPARE_MULTI}" \
  --run "${RUN1_LABEL}=${RUN1_OUT}" \
  --run "${RUN2_LABEL}=${RUN2_OUT}" \
  --run "${RUN3_LABEL}=${RUN3_OUT}" \
  --outdir "${COMPARE_OUT}"

echo "[INFO] =========================================================="
echo "[INFO] Merging finalised run"
echo "[INFO] =========================================================="

python "${MERGE_SAFE}" \
  --base-run "${RUN1_OUT}" \
  --add-run "${RUN2_LABEL}=${RUN2_OUT}" \
  --add-run "${RUN3_LABEL}=${RUN3_OUT}" \
  --outdir "${MERGED_OUT}"

echo "[INFO] =========================================================="
echo "[INFO] Summarising prmature stops from finalised runs"
echo "[INFO] =========================================================="

python "${PREMSTOP_SUMMARY}" \
  --run "${RUN1_LABEL}=${RUN1_OUT}" \
  --run "${RUN2_LABEL}=${RUN2_OUT}" \
  --run "${RUN3_LABEL}=${RUN3_OUT}" \
  --outdir "${PREMSTOP_OUT}"

echo "[INFO] Done"

echo "[INFO] Final run directories:"
echo "[INFO]   ${RUN1_OUT}"
echo "[INFO]   ${RUN2_OUT}"
echo "[INFO]   ${RUN3_OUT}"
echo "[INFO] Comparison output: ${COMPARE_OUT}"
echo "[INFO] Merged output:     ${MERGED_OUT}"
echo "[INFO] Premstop summary:  ${PREMSTOP_OUT}"