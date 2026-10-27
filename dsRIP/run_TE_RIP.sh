#!/bin/bash -l
#PBS -q cpu-test
#PBS -l nodes=1:ppn=8
#PBS -l mem=30gb
#PBS -N TEtranscripts_run
#PBS -o TE_out.txt
#PBS -e TE_err.txt

cd $PBS_O_WORKDIR

module purge
module load bioinfo/gcc-13.1.0
module load libs/gcc/zlib-1.2.13

source ~/.bashrc
conda activate TEtrans
unset PERL5LIB

TEtranscripts \
  -t \
    /iitjfs/home/2024rbb1003/projects/ASE_NB_project/RIPSeq/aligned_TE/BK846-01L0014/BK846-01L0014_Aligned.sortedByCoord.out.bam \
    /iitjfs/home/2024rbb1003/projects/ASE_NB_project/RIPSeq/aligned_TE/BK846-01L0015/BK846-01L0015_Aligned.sortedByCoord.out.bam \
    /iitjfs/home/2024rbb1003/projects/ASE_NB_project/RIPSeq/aligned_TE/BK846-01L0017/BK846-01L0017_Aligned.sortedByCoord.out.bam \
    /iitjfs/home/2024rbb1003/projects/ASE_NB_project/RIPSeq/aligned_TE/BK846-01L0018/BK846-01L0018_Aligned.sortedByCoord.out.bam \
  -c \
    /iitjfs/home/2024rbb1003/projects/ASE_NB_project/RIPSeq/aligned_TE/BK846-01L0013/BK846-01L0013_Aligned.sortedByCoord.out.bam \
    /iitjfs/home/2024rbb1003/projects/ASE_NB_project/RIPSeq/aligned_TE/BK846-01L0016/BK846-01L0016_Aligned.sortedByCoord.out.bam \
  --GTF /iitjfs/home/iitjmu11148/dbs/hg38/release_113/ensembl_annotation/Homo_sapiens.GRCh38.113.chr.gtf \
  --TE /home/2024rbb1003/projects/ASE_NB_project/RIPSeq/GRCh38_Ensembl_rmsk_TE.gtf \
  --sortByPos \
  --mode multi \
  --project RIPSeq_TE \
  --outdir TE_output_input

