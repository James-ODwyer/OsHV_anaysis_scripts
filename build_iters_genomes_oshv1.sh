#!/bin/bash
#SBATCH --account=XXXXX
#SBATCH --job-name=refine_consensus
#SBATCH --cpus-per-task=40
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 1
#SBATCH --mem=120G
#SBATCH --time=32:00:00


module load bowtie2/2.4.5
module load samtools/1.16.1
module load parallel



# Removed pathings from full script (put your paths here to replicate)

RUNDIR="/scratch3/*/OsHVrun"
GENOME_DIR="{RUNDIR}/genomes_rd1"
READ_DIR="{RUNDIR}/fastp"
OUT_DIR="{RUNDIR}/genomes_rd2"

mkdir -p $OUT_DIR


process_sample() {

  GENOME=$1
  BASE=$(basename $GENOME _genome_draft1.fa)

  R1=${READ_DIR}/${BASE}_filtered_R1.fastq.gz
  R2=${READ_DIR}/${BASE}_filtered_R2.fastq.gz

  SAMPLE_DIR=${OUT_DIR}/${BASE}
  mkdir -p $SAMPLE_DIR

  echo "Processing $BASE"

  bowtie2-build $GENOME ${SAMPLE_DIR}/${BASE}_index

  #iter 1

  bowtie2 \
    --very-sensitive \
    -x ${SAMPLE_DIR}/${BASE}_index \
    -1 $R1 \
    -2 $R2 \
    -S ${SAMPLE_DIR}/iter1.sam \
    -p 8

  samtools view -bS ${SAMPLE_DIR}/iter1.sam | \
    samtools sort -o ${SAMPLE_DIR}/iter1_sorted.bam

  samtools index ${SAMPLE_DIR}/iter1_sorted.bam

  # Extract mapped reads (-F 4 = mapped)
  samtools view -b -F 4 ${SAMPLE_DIR}/iter1_sorted.bam > ${SAMPLE_DIR}/iter1_mapped.bam

  # con 1

  samtools consensus \
    -d 1 \
    --call-fract 0.5 \
    ${SAMPLE_DIR}/iter1_mapped.bam > ${SAMPLE_DIR}/consensus1.fa


  bowtie2-build ${SAMPLE_DIR}/consensus1.fa ${SAMPLE_DIR}/consensus1_index

  #iter 2

  bowtie2 \
    --very-sensitive \
    -x ${SAMPLE_DIR}/consensus1_index \
    -1 $R1 \
    -2 $R2 \
    -S ${SAMPLE_DIR}/iter2.sam \
    -p 8

  samtools view -bS ${SAMPLE_DIR}/iter2.sam | \
    samtools sort -o ${SAMPLE_DIR}/iter2_sorted.bam

  samtools index ${SAMPLE_DIR}/iter2_sorted.bam

  samtools view -b -F 4 ${SAMPLE_DIR}/iter2_sorted.bam > ${SAMPLE_DIR}/iter2_mapped.bam


  # con 2

  samtools consensus \
    -d 5 \
    --call-fract 0.66 \
    ${SAMPLE_DIR}/iter2_mapped.bam > ${SAMPLE_DIR}/consensus2_final.fa

  echo "Finished $BASE"
}

export -f process_sample
export GENOME_DIR READ_DIR OUT_DIR

#Run 5 parallel samples, bowtie + smatools are very memory inexpensive 
find $GENOME_DIR -name "*_genome_draft1.fa" | \
  parallel -j 5 process_sample {}

echo "All samples complete"