eval "$(conda shell.bash hook)"
conda activate /home/vidh72t/.conda/envs/oshv1_annot



python /group/sequencing/assembly/James/OsHV_annot/prepare_genbank_feature_tables.py \
  --gff-dir /group/sequencing/assembly/James/OsHV_annot/Final_gffs \
  --fasta-dir /group/sequencing/assembly/James/OsHV_annot/oshv1_annotation_MF509813_primary/00_split_genomes \
  --outdir /group/sequencing/assembly/James/OsHV_annot/genbank_feature_tables \
  --gap-min 10 \
  --organism "Ostreid herpesvirus 1" \
  --mol-type "genomic DNA" \
  --taxon 26193 \
  --host "Pacific oyster" \
  --gap-type "within scaffold" \
  --linkage-evidence "paired-ends"
