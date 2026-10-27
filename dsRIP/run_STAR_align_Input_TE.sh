#!/bin/bash


sname=$1;
threads=$2;
module load bioinfo/miniconda3


inpath="/iitjfs/home/2024rbb1003/projects/ASE_NB_project/RIPSeq/Input/Data/";

outpath="/iitjfs/home/2024rbb1003/projects/ASE_NB_project/RIPSeq/aligned_TE/";

mkdir $outpath/${sname}

STAR \
  --runThreadN $threads \
  --readFilesCommand zcat \
  --genomeDir /iitjfs/home/iitjmu11148/dbs/hg38/release_113/STAR_ensembl_index151RL/ \
  --readFilesIn $inpath/Unknown_${sname}_good_1.fq.gz $inpath/Unknown_${sname}_good_2.fq.gz \
  --outFileNamePrefix $outpath/${sname}/${sname}_ \
  --outSAMtype BAM SortedByCoordinate \
  --quantMode GeneCounts \
  \
  --outFilterMultimapNmax 100 \
  --winAnchorMultimapNmax 100 \
  --outSAMmultNmax -1 \
  --outSAMattributes NH HI AS nM MD \
  --outFilterMismatchNmax 999 \
  --outFilterMismatchNoverLmax 0.04 \
  --alignEndsType EndToEnd ;


