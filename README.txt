OsHV-1 ANALYSIS SCRIPTS
=======================

Scripts supporting whole-genome characterisation of Australian and New Zealand
Ostreid herpesvirus 1 (OsHV-1) variants.


 Inputs and output locations are now supplied as command-line
"arguments. HPC module loading, Conda activation, SLURM account names, partitions,"
and email addresses are intentionally not embedded. Activate the appropriate
"software environment before running a script, or add site-specific scheduler"
headers to a local copy.

Repository layout
-----------------
CDS_annotation/
Data_filtering_and_assembly/
GARD/
Mapping_raws/
Sensitivity_test/
Similarity_analysis/

The six sections below correspond to these workflow areas.

General requirements
--------------------
"Core command-line tools used across the repository include Python 3, R, MAFFT,"
"IQ-TREE 2, HyPhy, fastp, Bowtie 2, SAMtools, SPAdes, BBMap, NanoFilt, Guppy,"
"Liftoff, minimap2, Prodigal, BLAST+, NCBI EDirect, and GNU-compatible shell"
utilities. Not every tool is required for every section.

"Python packages include Biopython, pandas, NumPy, SciPy, and matplotlib."
"R packages used to plot sensitivity trees include ape, phangorn, ggplot2,"
"ggtree, dplyr, stringr, tibble, and optionally svglite."

Recommended practice:
  1. Create a dedicated Conda/Mamba environment.
  2. Confirm each executable is available in PATH.
  3. Run scripts
  4. Use a new output directory for each analysis.


File naming assumptions are documented below. Sample identifiers should be
"consistent across FASTA, FASTQ, GFF3, and Newick files."


1. CDS_annotation
=================

Purpose
-------
"Reference-guided and homology-assisted annotation of OsHV-1 genomes, comparison"
"of annotations generated from multiple primary references, correction of gene"
"models affected by ambiguous sequence, recovery of plausible partial genes,"
"manual application of curated boundary decisions, summary of premature stop"
"codons, and preparation of NCBI five-column feature tables."

Main files
----------
oshv1_pipeline.sh
  Downloads reference GenBank records or annotates a multi-FASTA input. The
"  annotation mode splits genomes, parses references, applies Liftoff, predicts"
"  ORFs with Prodigal, uses protein homology to recover missed CDSs, and writes"
  per-sample final GFF3 files.

Run_oshv1_annotation_pipeline_v2.sh
"  Orchestrates three primary-reference annotation runs, truncate-aware"
"  correction, partial-gene recovery, cross-run comparison, safe merging, and"
  premature-stop summaries. The first argument is a project output directory.
  Optional second and third arguments override the input genome multi-FASTA and
  reference GBFF directory.

truncate_aware_gff_v2.py
"  Trims terminal ambiguous sequence, truncates CDSs at the first internal stop,"
"  marks partial boundaries, and optionally flags/removes likely start-gap"
  artefacts using reference CDS lengths.

recover_gap_partial_genes_v2.py
  Removes gap-associated short annotations and attempts to recover partial gene
  models by chained TBLASTN protein homology.

compare_oshv1_multi_runs.py
  Compares gene presence among two or more primary-reference runs and writes
"  per-reference, pairwise, multi-run, and summary tables."

merge_multiref_gffs_safe.py
  Uses one annotation run as the base and adds genuinely missing loci from
  additional runs while avoiding replacement of existing truncated models or
  import of overlapping duplicate annotations.

summarize_premature_stops_v2.py
  Summarises premature-stop candidates across multiple annotation runs while
  suppressing a candidate if another reference-based run supports a complete
  CDS for the same sample and gene.

apply_manual_gff_redo.py
  Applies manually curated per-sample and per-gene GFF3 replacements held in a
  CSV decision table.

update_bad_annotations.sh
  Portable wrapper for apply_manual_gff_redo.py.

prepare_genbank_feature_tables.py
"  Converts curated GFF3 and FASTA files to NCBI five-column feature tables,"
  adding source and assembly-gap records and preserving partial CDS notation. Used downstream for preparing tables for NCBI upload

create_feature_tables.sh
  Portable wrapper supplying the OsHV-1 defaults used by
  prepare_genbank_feature_tables.py. Used downstream for preparing tables for NCBI upload

Example commands
----------------
Download reference records:

  ./CDS_annotation/oshv1_pipeline.sh download_refs \
    -o data/reference_gbff

Run a single primary-reference annotation:

  ./CDS_annotation/oshv1_pipeline.sh annotate \
    -i data/oshv1_genomes.fasta \
    -r data/reference_gbff \
    -p MF509813.1 \
    -o results/annotation_MF509813 \
    -t 16 \
    --min-ident 75 \
    --min-qcov 75

Run the three-references separately:

  ./CDS_annotation/Run_oshv1_annotation_pipeline_v2.sh \
    results/annotation_project \
    data/oshv1_genomes.fasta \
    data/reference_gbff

Compare pre-existing annotation runs:

  python3 CDS_annotation/compare_oshv1_multi_runs.py \
    --run MF509813=results/annotation_MF509813 \
    --run NC_005881=results/annotation_NC_005881 \
    --run KY242785=results/annotation_KY242785 \
    --outdir results/annotation_comparison

Apply a manual GFF correction table:
* Note manual checking was performed before this step using Geneious Prime
  ./CDS_annotation/update_bad_annotations.sh \
    results/merged_annotations \
    data/manual_annotation_decisions.csv \
    results/merged_annotations_curated

Create GenBank feature tables:

  ./CDS_annotation/create_feature_tables.sh \
    --gff-dir results/final_gffs \
    --fasta-dir results/split_genomes \
    --outdir results/genbank_feature_tables

Important outputs
-----------------
"Per-sample final GFF3 files, annotation manifests, missing-gene tables,"
"pairwise run comparisons, merged GFF3 files, partial-gene recovery logs,"
"premature-stop summaries, and NCBI .tbl files."


2. Data_filtering_and_assembly
==============================

Purpose
-------
"Quality filtering of paired-end Illumina reads, host-read subtraction, de novo"
"SPAdes assembly, iterative read remapping and consensus refinement, targeted"
"Nanopore basecalling/filtering, and whole-genome phylogenetic inference."

Main files
----------
Run_prelim_host_depletion_and_assembly_OsHV.sh
"  Runs fastp, Bowtie 2 host subtraction, SAMtools extraction of unmapped read"
"  pairs, and SPAdes --careful assembly. Input reads must be named"
  SAMPLE_R1.fastq.gz and SAMPLE_R2.fastq.gz. The host Bowtie 2 index must be
  constructed before this script is run.

build_iters_genomes_oshv1.sh
  Maps filtered Illumina reads to draft genomes and generates two simple-mode
"  SAMtools consensuses. The depth-1, call-fraction-0.50 consensus is an"
  intermediate remapping sequence. The final consensus uses depth 5 and a
  call fraction of 0.66. Drafts are expected as SAMPLE_genome_draft1.fa and
  reads as SAMPLE_filtered_R1.fastq.gz/SAMPLE_filtered_R2.fastq.gz.

MinION_analysis.sh
"  Calls and trims Nanopore reads with a user-supplied Guppy executable,"
"  concatenates barcode FASTQs, and filters with NanoFilt. Defaults reproduce"
"  the analysis thresholds: Q12, minimum length 800 bp, maximum length 5000 bp."

MAFFt_IQtree.sh
  Aligns whole genomes with MAFFT using the manuscript parameters and runs an
  IQ-TREE ModelFinder analysis with SH-aLRT and ultrafast bootstrap support.

Example commands
----------------
Short-read preprocessing and assembly:

  ./Data_filtering_and_assembly/Run_prelim_host_depletion_and_assembly_OsHV.sh \
    --reads-dir data/illumina \
    --host-index data/host_index/pacific_oyster \
    --outdir results/preliminary_assembly \
    --threads 8 \
    --jobs 3

Iterative consensus refinement:

  ./Data_filtering_and_assembly/build_iters_genomes_oshv1.sh \
    --genome-dir data/draft_genomes \
    --read-dir results/preliminary_assembly/fastp \
    --outdir results/refined_genomes \
    --threads 8 \
    --jobs 5

Targeted Nanopore processing:

  ./Data_filtering_and_assembly/MinION_analysis.sh \
    --input-dir data/fast5 \
    --outdir results/minion \
    --guppy-bin /path/to/guppy_basecaller \
    --device auto

Whole-genome alignment and tree:

  ./Data_filtering_and_assembly/MAFFt_IQtree.sh \
    data/whole_genomes.fasta \
    results/whole_genome_phylogeny \
8


3. GARD
=======

Purpose
-------
"Screen a whole-genome alignment for putative changes in phylogenetic history,"
"split the alignment at GARD breakpoints, infer partition-specific trees, match"
"taxon sets before comparison, and quantify concordance using conventional and"
complementary tree metrics.

Main files
----------
Run_oshv1_GARD.sh
  Aligns genomes with MAFFT and launches HyPhy GARD.



gard_partition_tree_compare_v2.py
"  Splits an existing alignment at supplied breakpoints, removes taxa exceeding"
"  the partition missingness threshold, runs IQ-TREE, midpoint-roots plotting"
"  copies, compares taxon-matched partition and baseline trees, and reports"
  unfiltered and support-filtered RF/split statistics.

augment_gard_tree_comparisons.py
"  Adds Spearman correlation of patristic distances, least-squares branch-length"
"  scaling, scaled Pearson correlation and NRMSE, nearest-neighbour retention,"
"  quartet concordance, per-partition scatter plots, and a multi-metric summary."

Run_split_trees.sh
  Portable wrapper that runs the partition comparison and extended metric
  analysis with the study defaults.

Example commands
----------------
Run GARD:

  ./GARD/Run_oshv1_GARD.sh \
    data/gard_input.fasta \
    results/gard/oshv1 \
    36 \
18

Compare the resulting partitions:
*Note these are the positions identified as breakpoints in the above study
  ./GARD/Run_split_trees.sh \
    results/gard/oshv1.aligned.fasta \
"    20253,46357,67872,89549,178578,182545,186115,188877,216228 \"
    results/gard/baseline.treefile \
    results/gard/partition_comparison

Interpretation guidance
-----------------------
RF distance is retained as a conventional exact-split metric. Supported-split
precision distinguishes strongly supported contradiction from reduced
resolution. High precision with lower recall indicates that a partition
supports a subset of the baseline relationships rather than many alternative
"relationships. Patristic-distance correlation, scaled branch-length agreement,"
and quartet concordance provide complementary assessments of similarity.


4. Mapping_raws
===============

Purpose
-------
"Map sample reads to complete genomes with BBMap, distribute equally scoring"
"ambiguous reads across repeated regions, extract mapped reads, calculate depth,"
and create individual and combined coverage plots.
Raw reads can be found on NCBI Sequence Read Archive.

Main files
----------
map_raws_to_genomes.sh
  Performs one-to-one sample-to-genome mapping. Genome filenames and read
  prefixes must share the same sample identifier as the reads and genomes are paired through regex pattern matching.


plot_coverage_oshv.py
  Reads BBMap pileup per-base tables and generates one coverage plot per sample
"  plus a multi-panel overview. Input and output directories, filename pattern,"
  and coverage threshold are command-line options.

Example commands
----------------
Map reads and generate depth tables:

  ./Mapping_raws/map_raws_to_genomes.sh \
    --genome-dir data/final_genomes \
    --reads-dir data/viral_reads \
    --outdir results/read_mapping \
    --threads 8 \
    --minid 0.8

Plot coverage:

  python3 Mapping_raws/plot_coverage_oshv.py \
    --input-dir results/read_mapping \
    --outdir results/read_mapping/coverage_plots \
    --threshold 5

Cross-sample mappings used for specific quality-control questions should be run
as separate explicit commands with clearly named output directories rather than
being hard-coded into the general one-to-one mapping script.


5. Sensitivity_test
===================

Purpose
-------
Construct a partitioned phylogeny from complete CDS regions meeting the
"minimum-depth-10 sensitivity criterion, plot the resulting tree, build a"
"matched whole-genome baseline, and quantify agreement between the whole-genome"
and high-coverage CDS analyses.

Main files
----------
align_and_run_partitioned_iqtree_dp10.sh
"  Aligns each retained CDS independently, concatenates alignments by taxon,"
"  inserts Ns for missing loci, writes a NEXUS partition file, and infers a"
  partitioned IQ-TREE phylogeny with a separate model selected per CDS.

maketrees.R
  Publication figure generator for the DP>=10 tree. Requires a tree path and
  output directory. An optional third argument is a comma-separated list of
  labels to omit from the plotted/exported tree.

run_maketrees.sh
  Portable wrapper for maketrees.R.

run_dp10_cds_vs_whole_genome_comparison.sh
  Creates a whole-genome alignment and baseline tree from the same taxa present
"  in the DP>=10 CDS tree, then reports RF, support-filtered split, patristic,"
"  nearest-neighbour, quartet, and scaled branch-length metrics and figures."

Example commands
----------------
Build the DP>=10 partitioned CDS tree:

  ./Sensitivity_test/align_and_run_partitioned_iqtree_dp10.sh \
    data/dp10_complete_cds_by_gene \
    results/dp10_partitioned_tree \
16

Plot the tree:

  ./Sensitivity_test/run_maketrees.sh \
    results/dp10_partitioned_tree/concatenated/OsHV_dp10_partitioned.treefile \
    results/dp10_partitioned_tree/figures

Compare with a taxon-matched whole-genome baseline:

  ./Sensitivity_test/run_dp10_cds_vs_whole_genome_comparison.sh \
    results/dp10_partitioned_tree \
    data/all_whole_genomes.fasta \
    results/dp10_partitioned_tree/whole_genome_comparison

"The key interpretation is whether the Australian and New Zealand grouping,"
"including the 2024 Tasmanian cluster, is retained after restricting the"
analysis to complete CDS sequence supported at DP>=10. Changes in the relative
position of external references may reflect the different genomic components
represented by whole-genome and stringently filtered coding datasets.


6. Similarity_analysis
======================

Purpose
-------
This folder is reserved for scripts that produce pairwise or sliding-window
whole-genome similarity profiles used to localise genomic divergence among
representative OsHV-1 genomes.

No similarity-analysis source file was included in the uploaded batch used to
prepare this sanitised release. Add the relevant script to this folder before
"publication, ensure that input genomes and output paths are command-line"
"arguments, and document the window size, step size, treatment of gaps, reference"
"selection, and coordinate system. The intended analysis should be described as"
"localising and quantifying regional similarity, not as formal structural-variant"
detection unless a dedicated structural-variant method is applied.

Suggested command-line interface:

  python3 Similarity_analysis/create_similarity_plots_OSHV.py \
    --alignment data/whole_genome_alignment.fasta \
    --focal-samples data/australian_new_zealand_samples.txt \
    --references data/reference_samples.txt \
    --window-size 1000 \
    --step-size 100 \
    --outdir results/similarity_profiles


Reproducibility and portability notes
=====================================

1. The  files do not activate a named Conda environment or load
   institution-specific modules. Environment setup belongs in a cluster job
   submission wrapper or documented environment file.
"2. Scheduler resource requests were removed. CPU, memory, walltime, account,"
"   partition, and notification settings should be supplied by the user."
3. Absolute project paths were replaced by positional or named command-line
   arguments.
4. Shell scripts use strict mode where practical and quote paths.
"5. Before analysing production data, test each workflow on a small subset and"
   inspect logs and intermediate files. The analyses scripts are not extensively tested for edge cases or data outside of the utilised study.
6. A conda environment yaml file is provided which contains the majority of packages used.
" however, some programs and packages were loaded as modules and so will need to either be "
loaded on your HPC system  or installed via conda into the working environment.
