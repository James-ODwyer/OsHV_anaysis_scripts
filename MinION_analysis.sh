#!/bin/bash
#SBATCH --account=XXXXX
#SBATCH --job-name=minion_basecall_filter
#SBATCH --cpus-per-task=12
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 1
#SBATCH --mem=40G
#SBATCH --time=48:00:00

module load NanoFilt/2.8.0

RUNDIR="/scratch3/*/OsHVrun"
RAW_DIR="/scratch3/*/OsHVraws_MinION"
OUT_DIR="${RUNDIR}/minion_processed"

GUPPY_OUT="${OUT_DIR}/basecalled"
FILTERED_OUT="${OUT_DIR}/filtered"

mkdir -p $GUPPY_OUT
mkdir -p $FILTERED_OUT


# Guppy locally downloaded onto working directory
/localdir/guppy_basecaller \
  --enable_trim_barcodes \
  --trim_adapters \
  --barcode_kits "SQK-RPB004" \
  -i ${RAW_DIR} \
  -s ${GUPPY_OUT} \
  -c dna_r9.4.1_450bps_hac.cfg \
  --device cuda:${CUDA_VISIBLE_DEVICES}

for BARCODE_DIR in ${GUPPY_OUT}/barcode*; do
  BARCODE=$(basename $BARCODE_DIR)

  cat ${BARCODE_DIR}/*.fastq > ${GUPPY_OUT}/${BARCODE}.fastq
done


for FASTQ in ${GUPPY_OUT}/barcode*.fastq; do
  BASE=$(basename $FASTQ .fastq)

  NanoFilt \
    -l 800 \
    --maxlength 5000 \
    < $FASTQ > ${FILTERED_OUT}/${BASE}_filtered.fastq
done
