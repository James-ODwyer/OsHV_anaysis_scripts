#!/bin/bash
#SBATCH --nodes=1
#SBATCH --job-name="bbmapany"
#SBATCH --account="aggstrategic1"
#SBATCH --partition="batch"
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=48GB
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=james.o\'dwyer@agriculture.vic.gov.au
#SBATCH --time=48:0:00





module load BBMap/
module load SAMtools/


reflist=(`ls /group/sequencing/assembly/James/OsHV_annot/Align_raws/genomes/*`) 
readsdir=/group/sequencing/assembly/James/OsHV_annot/Align_raws 


outdir=/group/sequencing/assembly/James/OsHV_annot/Align_raws/bbmapaligned_to_final_genomes_for_SRA_upload_1to1_2


mkdir "$outdir"



for file in ${reflist[@]}; do

name=(`awk '{sub(/.fsa.*/, ""); sub(/.*genomes\//, ""); print }'<<< $file`)


R1="$readsdir"/"$name"_R1_from_hits.fastq.gz
R2="$readsdir"/"$name"_R2_from_hits.fastq.gz

bbmap.sh in1="$R1" in2="$R2" ref="$file" out="$outdir"/"$name"_aligned.sam minid=0.8 threads=8 ambig=random
samtools view -h -F 4 "$outdir"/"$name"_aligned.sam > "$outdir"/"$name"_aligned_hits.sam


	# Double check if fails #
  samtools fastq "$outdir"/"$name"_aligned_hits.sam \
    -1 "$outdir"/"$name"_R1_from_hits.fastq.gz \
    -2 "$outdir"/"$name"_R2_from_hits.fastq.gz \
    -0 /dev/null -s /dev/null -n
#




pileup.sh in="$outdir"/"$name"_aligned.sam out="$outdir"/"$name"_depth_stats.txt basecov="$outdir"/"$name"_perbase_stats

done




outdir2=/group/sequencing/assembly/James/OsHV_annot/Align_raws/bbmapaligned_to_final_genomes_for_SRA_upload_mappingtohighcov
mkdir "$outdir2"


name="OsHV_NSW_2010_1"
reffile="/group/sequencing/assembly/James/OsHV_annot/Align_raws/genomes/OsHV_NSW_2011_1.fsa"
R1="$readsdir"/"$name"_R1_from_hits.fastq.gz
R2="$readsdir"/"$name"_R2_from_hits.fastq.gz

bbmap.sh in1="$R1" in2="$R2" ref="$reffile" out="$outdir2"/"$name"_aligned.sam minid=0.8 threads=8 ambig=random
samtools view -h -F 4 "$outdir2"/"$name"_aligned.sam > "$outdir2"/"$name"_aligned_hits.sam


	# Double check if fails #
  samtools fastq "$outdir2"/"$name"_aligned_hits.sam \
    -1 "$outdir2"/"$name"_R1_from_hits.fastq.gz \
    -2 "$outdir2"/"$name"_R2_from_hits.fastq.gz \
    -0 /dev/null -s /dev/null -n
#




pileup.sh in="$outdir2"/"$name"_aligned.sam out="$outdir2"/"$name"_depth_stats.txt basecov="$outdir2"/"$name"_perbase_stats




outdir2=/group/sequencing/assembly/James/OsHV_annot/Align_raws/bbmapaligned_to_final_genomes_for_SRA_upload_mappingtohighcov

name="OsHV_TAS_2024_2"
reffile="/group/sequencing/assembly/James/OsHV_annot/Align_raws/genomes/OsHV_TAS_2024_1.fsa"
R1="$readsdir"/"$name"_R1_from_hits.fastq.gz
R2="$readsdir"/"$name"_R2_from_hits.fastq.gz

bbmap.sh in1="$R1" in2="$R2" ref="$reffile" out="$outdir2"/"$name"_aligned.sam minid=0.8 threads=8 ambig=random
samtools view -h -F 4 "$outdir2"/"$name"_aligned.sam > "$outdir2"/"$name"_aligned_hits.sam


	# Double check if fails #
  samtools fastq "$outdir2"/"$name"_aligned_hits.sam \
    -1 "$outdir2"/"$name"_R1_from_hits.fastq.gz \
    -2 "$outdir2"/"$name"_R2_from_hits.fastq.gz \
    -0 /dev/null -s /dev/null -n
#




pileup.sh in="$outdir2"/"$name"_aligned.sam out="$outdir2"/"$name"_depth_stats.txt basecov="$outdir2"/"$name"_perbase_stats






outdir2=/group/sequencing/assembly/James/OsHV_annot/Align_raws/bbmapaligned_to_final_genomes_for_SRA_upload_mappingtohighcov

name="OsHV_TAS_2024_3"
reffile="/group/sequencing/assembly/James/OsHV_annot/Align_raws/genomes/OsHV_TAS_2024_1.fsa"
R1="$readsdir"/"$name"_R1_from_hits.fastq.gz
R2="$readsdir"/"$name"_R2_from_hits.fastq.gz

bbmap.sh in1="$R1" in2="$R2" ref="$reffile" out="$outdir2"/"$name"_aligned.sam minid=0.8 threads=8 ambig=random
samtools view -h -F 4 "$outdir2"/"$name"_aligned.sam > "$outdir2"/"$name"_aligned_hits.sam


	# Double check if fails #
  samtools fastq "$outdir2"/"$name"_aligned_hits.sam \
    -1 "$outdir2"/"$name"_R1_from_hits.fastq.gz \
    -2 "$outdir2"/"$name"_R2_from_hits.fastq.gz \
    -0 /dev/null -s /dev/null -n
#




pileup.sh in="$outdir2"/"$name"_aligned.sam out="$outdir2"/"$name"_depth_stats.txt basecov="$outdir2"/"$name"_perbase_stats




outdir3=/group/sequencing/assembly/James/OsHV_annot/Align_raws/bbmapaligned_tas2024_against_2011NSW

mkdir "$outdir3"

name="OsHV_TAS_2024_1"
reffile="/group/sequencing/assembly/James/OsHV_annot/Align_raws/genomes/OsHV_NSW_2011_1.fsa"
R1="$readsdir"/"$name"_R1_from_hits.fastq.gz
R2="$readsdir"/"$name"_R2_from_hits.fastq.gz

bbmap.sh in1="$R1" in2="$R2" ref="$reffile" out="$outdir3"/"$name"_aligned.sam minid=0.8 threads=8 ambig=random
samtools view -h -F 4 "$outdir3"/"$name"_aligned.sam > "$outdir3"/"$name"_aligned_hits.sam


	# Double check if fails #
  samtools fastq "$outdir3"/"$name"_aligned_hits.sam \
    -1 "$outdir3"/"$name"_R1_from_hits.fastq.gz \
    -2 "$outdir3"/"$name"_R2_from_hits.fastq.gz \
    -0 /dev/null -s /dev/null -n
#




pileup.sh in="$outdir3"/"$name"_aligned.sam out="$outdir3"/"$name"_depth_stats.txt basecov="$outdir3"/"$name"_perbase_stats



outdir3=/group/sequencing/assembly/James/OsHV_annot/Align_raws/bbmapaligned_tas2024_against_2011NSW

name="OsHV_TAS_2024_2"
reffile="/group/sequencing/assembly/James/OsHV_annot/Align_raws/genomes/OsHV_NSW_2011_1.fsa"
R1="$readsdir"/"$name"_R1_from_hits.fastq.gz
R2="$readsdir"/"$name"_R2_from_hits.fastq.gz

bbmap.sh in1="$R1" in2="$R2" ref="$reffile" out="$outdir3"/"$name"_aligned.sam minid=0.8 threads=8
samtools view -h -F 4 "$outdir3"/"$name"_aligned.sam > "$outdir3"/"$name"_aligned_hits.sam


	# Double check if fails #
  samtools fastq "$outdir3"/"$name"_aligned_hits.sam \
    -1 "$outdir3"/"$name"_R1_from_hits.fastq.gz \
    -2 "$outdir3"/"$name"_R2_from_hits.fastq.gz \
    -0 /dev/null -s /dev/null -n
#




pileup.sh in="$outdir3"/"$name"_aligned.sam out="$outdir3"/"$name"_depth_stats.txt basecov="$outdir3"/"$name"_perbase_stats






outdir3=/group/sequencing/assembly/James/OsHV_annot/Align_raws/bbmapaligned_tas2024_against_2011NSW

name="OsHV_TAS_2024_3"
reffile="/group/sequencing/assembly/James/OsHV_annot/Align_raws/genomes/OsHV_NSW_2011_1.fsa"
R1="$readsdir"/"$name"_R1_from_hits.fastq.gz
R2="$readsdir"/"$name"_R2_from_hits.fastq.gz

bbmap.sh in1="$R1" in2="$R2" ref="$reffile" out="$outdir3"/"$name"_aligned.sam minid=0.8 threads=8 ambig=random
samtools view -h -F 4 "$outdir3"/"$name"_aligned.sam > "$outdir3"/"$name"_aligned_hits.sam


	# Double check if fails #
  samtools fastq "$outdir3"/"$name"_aligned_hits.sam \
    -1 "$outdir3"/"$name"_R1_from_hits.fastq.gz \
    -2 "$outdir3"/"$name"_R2_from_hits.fastq.gz \
    -0 /dev/null -s /dev/null -n
#



pileup.sh in="$outdir3"/"$name"_aligned.sam out="$outdir3"/"$name"_depth_stats.txt basecov="$outdir3"/"$name"_perbase_stats



