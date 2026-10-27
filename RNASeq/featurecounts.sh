#!/bin/bash
#PBS -q cpu-test
#PBS -N Featurecounts
#PBS -l nodes=1:ppn=20
#PBS -l mem=10gb
#PBS -o xFeaturecounts_out.txt
#PBS -e xFeaturecounts_err.txt

module load bioinfo/subread-2.0.2

#featureCounts \
#-p -T 6 -t exon -g gene_id -Q 30 -F GTF \
#-a /iitjfs/home/iitjmu11148/dbs/hg38/release_113/ensembl_annotation/Homo_sapiens.GRCh38.113.chr.gtf \
#-o /iitjfs/home/2024rbb1003/projects/ASE_NB_project/RNASeq/Quantification_gene/hg38_counts.txt \
#$(cat /iitjfs/home/2024rbb1003/projects/ASE_NB_project/RNASeq/scripts/bam_list.txt)

featureCounts \
-p -T 6 -t exon -g transcript_id -Q 30 -F GTF \
-a /iitjfs/home/iitjmu11148/dbs/hg38/release_113/ensembl_annotation/Homo_sapiens.GRCh38.113.chr.gtf \
-o /iitjfs/home/2024rbb1003/projects/ASE_NB_project/RNASeq/Quantification_gene/hg38_counts_transcript_level.txt \
$(cat /iitjfs/home/2024rbb1003/projects/ASE_NB_project/RNASeq/scripts/bam_list.txt)
