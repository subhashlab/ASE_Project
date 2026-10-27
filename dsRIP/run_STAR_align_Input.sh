#!/bin/bash


sname=$1;
threads=$2;
module load bioinfo/miniconda3


inpath="/iitjfs/home/iitjmu11148/projects/Meena_project/ASE_proj/RIPseq/Input/Data";

outpath="/iitjfs/home/2024rbb1003/projects/ASE_NB_project/RIPSeq/aligned/";

mkdir $outpath/${sname}

STAR --runThreadN $threads --readFilesCommand zcat --genomeDir /iitjfs/home/iitjmu11148/dbs/hg38/release_113/STAR_ensembl_index151RL/ --readFilesIn $inpath/Unknown_${sname}_good_1.fq.gz $inpath/Unknown_${sname}_good_2.fq.gz --outFileNamePrefix $outpath/${sname}/${sname} --outSAMtype BAM SortedByCoordinate --quantMode GeneCounts ;


