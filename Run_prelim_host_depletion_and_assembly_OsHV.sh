#!/bin/bash 
#SBATCH --account=OD-229285
#SBATCH --job-name OsHV_prefilter_genome_builds
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 1
#SBATCH --cpus-per-task 40
#SBATCH --mem 126G
#SBATCH --time 108:00:00


#modules
module load fastp/0.23.2
module load bowtie2/2.4.5
module load samtools/1.16.1 
module load spades/3.15.5
module load parallel


# Removed pathings from full script (put your paths here to replicate)
READ_DIR="/scratch3/*/OsHVraws"


RUNDIR="/scratch3/*/OsHVrun"

FASTP_DIR="${RUNDIR}/fastp"
ALIGN_DIR="${RUNDIR}/bowtie2"
UNMAP_DIR="${RUNDIR}/unmapped"
SPADES_DIR="${RUNDIR}/spades"

REF="${RUNDIR}/GCA_025765685.3"

mkdir -p $FASTP_DIR $ALIGN_DIR $UNMAP_DIR $SPADES_DIR

# Create function for parallelisation across samples

process_sample() {
  R1=$1
  BASE=$(basename $R1 _R1.fastq.gz)
  R2=${READ_DIR}/${BASE}_R2.fastq.gz

  echo "Processing $BASE"

  # fastp
  fastp \
    -i $R1 \
    -I $R2 \
    -o ${FASTP_DIR}/${BASE}_filtered_R1.fastq.gz \
    -O ${FASTP_DIR}/${BASE}_filtered_R2.fastq.gz \
    --length_required 75 \
    --qualified_quality_phred 15 \
    --cut_front \
    --cut_front_window_size 4 \
    --cut_front_mean_quality 20 \
    --thread 8 \
    --html ${FASTP_DIR}/${BASE}_fastp.html \
    --json ${FASTP_DIR}/${BASE}_fastp.json

  # Bowtie2
  bowtie2 \
    --very-sensitive \
    -x $REF \
    -1 ${FASTP_DIR}/${BASE}_filtered_R1.fastq.gz \
    -2 ${FASTP_DIR}/${BASE}_filtered_R2.fastq.gz \
    -S ${ALIGN_DIR}/${BASE}.sam \
    -p 8

  samtools view -bS ${ALIGN_DIR}/${BASE}.sam > ${ALIGN_DIR}/${BASE}.bam
  samtools sort ${ALIGN_DIR}/${BASE}.bam -o ${ALIGN_DIR}/${BASE}_sorted.bam
  samtools index ${ALIGN_DIR}/${BASE}_sorted.bam

  # Collect unmapped as putative viral reads
  samtools view -b -f 12 -F 256 \
    ${ALIGN_DIR}/${BASE}_sorted.bam > ${UNMAP_DIR}/${BASE}_unmapped.bam

  samtools fastq \
    -1 ${UNMAP_DIR}/${BASE}_R1.fastq.gz \
    -2 ${UNMAP_DIR}/${BASE}_R2.fastq.gz \
    -0 /dev/null -s /dev/null -n \
    ${UNMAP_DIR}/${BASE}_unmapped.bam

  # Run spades, run as standard spades as there is almost exclusively viral reads in here now. 
  spades.py \
    -1 ${UNMAP_DIR}/${BASE}_R1.fastq.gz \
    -2 ${UNMAP_DIR}/${BASE}_R2.fastq.gz \
    --careful \
    -t 8 \
    -m 32 \
    -o ${SPADES_DIR}/${BASE}

  echo "Finished $BASE"
}

export -f process_sample
export READ_DIR FASTP_DIR ALIGN_DIR UNMAP_DIR SPADES_DIR REF

# 3 samples at a time as Spades can be memory intensuve

find $READ_DIR -name "*_R1.fastq.gz" | \
  parallel -j 3 process_sample {}

echo "Finished run"