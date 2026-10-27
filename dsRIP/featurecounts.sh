#!/bin/bash
#PBS -q cpu-test
#PBS -N ASE_QC
#PBS -l nodes=1:ppn=20
#PBS -l mem=10gb
#PBS -o xASE_out.txt
#PBS -e xASE_err.txt

module load bioinfo/subread-2.0.2

echo "#### QC started `date` \n";

featureCounts \
-p -T 6 -t exon -g transcript_id -Q 30 -F GTF \
-a /iitjfs/home/iitjmu11148/dbs/hg38/release_113/ensembl_annotation/Homo_sapiens.GRCh38.113.chr.gtf \
-o /iitjfs/home/2024rbb1003/projects/ASE_NB_project/RIPSeq/quantification/hg38_counts_transcript_level.txt \
$(cat /iitjfs/home/2024rbb1003/projects/ASE_NB_project/RIPSeq/quantification/bam_list.txt)

echo "#### QC completed `date` \n";
wait
