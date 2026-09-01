#!/bin/bash
#SBATCH --nodes=1
#SBATCH --job-name="GARD_genomes_trees"
#SBATCH --account="aggstrategic1"
#SBATCH --partition="batch"
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=70GB
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=james.o\'dwyer@agriculture.vic.gov.au
#SBATCH --time=120:0:00


#module load IQ-TREE



eval "$(conda shell.bash hook)"
conda activate /home/vidh72t/.conda/envs/oshv1_annot





#iqtree \
#  -s /group/sequencing/assembly/James/OsHV_annot/GARD/rebuild_alignmentfirst/aligned.fsa \
#  -b 2000 \
#  -alrt 2000 \
#  -m MFP \
#  -T AUTO \
#  --threads-max 24





python3 gard_partition_tree_compare_v2.py \
  --alignment /group/sequencing/assembly/James/OsHV_annot/GARD/rebuild_alignmentfirst/aligned.fsa \
  --breakpoints 20253,46357,67872,89549,178578,182545,186115,188877,216228 \
  --baseline-tree /group/sequencing/assembly/James/OsHV_annot/GARD/rebuild_alignmentfirst/IQtree_first/aligned.fsa.treefile \
  --outdir /group/sequencing/assembly/James/OsHV_annot/GARD/tree_compare_10bp_midpoint_20260702_v3 \
  --run-iqtree \
  --threads 24 \
  --iqtree-exe iqtree \
  --iqtree-model MFP \
  --ufboot 1000 \
  --iqtree-extra "--alrt 1000" \
  --drop-missing-above 0.5 \
  --min-taxa 4 \
  --min-usable-sites 20 \
  --support-thresholds 0.5,0.8,0.9 \
  --paired-support-mode min


python3 augment_gard_tree_comparisons.py \
  --summary /group/sequencing/assembly/James/OsHV_annot/GARD/tree_compare_10bp_midpoint_20260702_v3/partition_tree_summary.tsv \
  --outdir /group/sequencing/assembly/James/OsHV_annot/GARD/tree_compare_10bp_midpoint_20260702_v3/extended_tree_metrics