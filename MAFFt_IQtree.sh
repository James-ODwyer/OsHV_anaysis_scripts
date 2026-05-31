#!/bin/bash
#SBATCH --account=XXXXX
#SBATCH --job-name=mafft_iqtree_genomes
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --time=48:00:00
#SBATCH --ntasks-per-node=1
#SBATCH --nodes=1

module load mafft/7.490
module load iqtree


RUNDIR="/path/to/run"
INPUT_FASTA="${RUNDIR}/whole_genomes.fasta"
ALIGN_DIR="${RUNDIR}/alignment"
TREE_DIR="${RUNDIR}/iqtree"

mkdir -p $ALIGN_DIR
mkdir -p $TREE_DIR

mafft \
  --auto \
  --maxiterate 1000 \
  --retree 150 \
  --op 2 \
  --ep 0.12 \
  --adjustdirection \
  --thread 8 \
  ${INPUT_FASTA} > ${ALIGN_DIR}/whole_genome_OSHV_alignment.mafft

iqtree \
  -s ${ALIGN_DIR}/whole_genome_OSHV_alignment.mafft \
  -b 2000 \
  -alrt 2000 \
  -T AUTO \
  --threads-max 8 \
  -o REFERENCE_NAME \
  -pre ${TREE_DIR}/genome_tree
