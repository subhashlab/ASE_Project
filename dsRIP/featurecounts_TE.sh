#!/bin/bash
#PBS -q cpu-test
#PBS -N ASE_QC
#PBS -l nodes=1:ppn=20
#PBS -l mem=10gb
#PBS -o xASE_out.txt
#PBS -e xASE_err.txt

module load bioinfo/subread-2.0.2

echo "#### QC started $(date)"

xargs -a /iitjfs/home/2024rbb1003/projects/ASE_NB_project/RIPSeq/quantification/bam_list_TE.txt \
featureCounts \
  -p \
  -T 6 \
  -t exon \
  -g transcript_id \
  -Q 30 \
  -F GTF \
  -M \
  -O \
  --fraction \
  -a /iitjfs/home/2024rbb1003/projects/ASE_NB_project/RIPSeq/GRCh38_Ensembl_rmsk_TE.gtf \
  -o /iitjfs/home/2024rbb1003/projects/ASE_NB_project/RIPSeq/quantification/hg38_counts_TE.txt

echo "#### QC completed $(date)"
wait

