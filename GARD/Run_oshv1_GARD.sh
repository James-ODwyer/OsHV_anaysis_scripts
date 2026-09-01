#!/bin/bash
#SBATCH --nodes=1
#SBATCH --job-name="GARD_genomes"
#SBATCH --account="acct"
#SBATCH --partition="batch"
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=36
#SBATCH --mem=90GB
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=user@user.com
#SBATCH --time=128:0:00

set -euo pipefail

eval "$(conda shell.bash hook)"
conda activate /*/.conda/envs/oshv1_annot

module load MAFFT


mafft \
  --auto \
  --maxiterate 1000 \
  --retree 150 \
  --op 2 \
  --ep 0.12 \
  --adjustdirection \
  --thread 16 \
  All_genomes_fastas_except_3_high_missing.fasta > aligned.fsa



hyphy CPU=36 GARD \
  --alignment aligned.fsa \
  --max-breakpoints 18 \
  --mode Faster \
  ENV="TOLERATE_NUMERICAL_ERRORS=1;"


