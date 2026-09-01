
module load R/4.4.2-gfbf-2024a
module load R-bundle-Bioconductor
module load R-bundle-CRAN



#Rscript maketrees.R \
#  /group/sequencing/assembly/James/OsHV_annot/Align_raws/dp10_complete_CDS_sets/mafft_iqtree_partitioned_dp10/concatenated/OsHV_dp10_partitioned.treefile \
#  /group/sequencing/assembly/James/OsHV_annot/Align_raws/dp10_complete_CDS_sets/mafft_iqtree_partitioned_dp10/figures
  
  
#Rscript maketrees.R \
#  /group/sequencing/assembly/James/OsHV_annot/Align_raws/dp10_complete_CDS_sets/mafft_iqtree_partitioned_dp10_added_SRA/concatenated/OsHV_dp10_partitioned.treefile \
#  /group/sequencing/assembly/James/OsHV_annot/Align_raws/dp10_complete_CDS_sets/mafft_iqtree_partitioned_dp10_added_SRA/figures

Rscript maketrees.R \
  /group/sequencing/assembly/James/OsHV_annot/Align_raws/dp10_complete_CDS_sets/mafft_iqtree_partitioned_dp10_added_SRA_nodups/concatenated/OsHV_dp10_partitioned.treefile \
  /group/sequencing/assembly/James/OsHV_annot/Align_raws/dp10_complete_CDS_sets/mafft_iqtree_partitioned_dp10_added_SRA_nodups/figures
