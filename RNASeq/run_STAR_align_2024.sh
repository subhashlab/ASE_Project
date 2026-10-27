#!/bin/bash


sname=$1;
threads=$2;
module load bioinfo/miniconda3


inpath="/iitjfs/home/iitjmu11148/projects/Meena_project/ASE_proj/FASTQ_2024/";

outpath="/iitjfs/home/2024rbb1003/projects/ASE_NB_project/RNASeq/aligned/aligned_2024";

mkdir $outpath/${sname}

STAR \
 --runThreadN $threads \
 --readFilesCommand zcat \
 --genomeDir /iitjfs/home/iitjmu11148/dbs/hg38/release_113/STAR_ensembl_index151RL/ \
 --readFilesIn \
   $inpath/${sname}_L1_1.fq.gz \
   $inpath/${sname}_L1_2.fq.gz \
 --outFileNamePrefix $outpath/${sname}/${sname}. \
 --outSAMtype BAM SortedByCoordinate \
 --outSAMstrandField intronMotif \
 --outFilterType BySJout \
 --alignSJoverhangMin 8 \
 --alignSJDBoverhangMin 1 \
 --outFilterMismatchNmax 999 \
 --outFilterMismatchNoverReadLmax 0.04 \
 --alignIntronMin 20 \
 --alignIntronMax 1000000 \
 --alignMatesGapMax 1000000 ;



