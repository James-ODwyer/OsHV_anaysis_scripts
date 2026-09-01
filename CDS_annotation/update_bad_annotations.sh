eval "$(conda shell.bash hook)"
conda activate /home/vidh72t/.conda/envs/oshv1_annot



python /group/sequencing/assembly/James/OsHV_annot/apply_manual_gff_redo.py \
  --merged-dir /group/sequencing/assembly/James/OsHV_annot/oshv1_annotation_merged_MF509813plus \
  --csv /group/sequencing/assembly/James/OsHV_annot/REdo_these_genes.csv \
  --outdir /group/sequencing/assembly/James/OsHV_annot/oshv1_annotation_merged_MF509813plus_manual
