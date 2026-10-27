#!/bin/bash
#PBS -q cpu-test
#PBS -N rMATS_cyto_trt_vs_ctrl
#PBS -l nodes=1:ppn=20
#PBS -l mem=10gb
#PBS -o rMATS_cyto_trt_vs_ctrl_out_2.txt
#PBS -e rMATS_cyto_trt_vs_ctrl_err_2.txt


echo "Job started on $(hostname) at $(date)"

unset PERL5LIB
unset PERL_LOCAL_LIB_ROOT
unset PERL_MB_OPT
unset PERL_MM_OPT

# Activate conda environment
source ~/.bashrc
conda activate rmats_env

# Run rMATS
rmats.py \
  --b1 /home/2024rbb1003/projects/ASE_NB_project/RNASeq/scripts/Cyto_trt.txt \
  --b2 /home/2024rbb1003/projects/ASE_NB_project/RNASeq/scripts/Cyto_ctrl.txt \
  --gtf /iitjfs/home/iitjmu11148/dbs/hg38/release_113/ensembl_annotation/Homo_sapiens.GRCh38.113.chr.gtf \
  -t paired \
  --readLength 150 \
  --nthread 16 \
  --od /home/2024rbb1003/projects/ASE_NB_project/RNASeq/rmats_cyto_trt_vs_ctrl \
  --tmp /home/2024rbb1003/projects/ASE_NB_project/RNASeq/rmats_cyto_trt_vs_ctrl/tmp \
  --libType fr-secondstrand \
  --allow-clipping

echo "Job finished at $(date)"
