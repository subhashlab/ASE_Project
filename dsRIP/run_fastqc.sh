#!/bin/bash
#PBS -q cpu-test
#PBS -N ASE_QC
#PBS -l nodes=1:ppn=20
#PBS -l mem=10gb
#PBS -o xASE_out.txt
#PBS -e xASE_err.txt

#inpath="/iitjfs/home/2024rbb1003/projects/ASE_NB_project/RIPSeq/Input/Data/";
inpath="/iitjfs/home/2024rbb1003/projects/ASE_NB_project/RIPSeq/RIP/BMK_DATA_20231110112317_1/";
#outpath="/iitjfs/home/2024rbb1003/projects/RIPSeq/QC/QC_Input";
outpath="/iitjfs/home/2024rbb1003/projects/ASE_NB_project/RIPSeq/QC/QC_RIP";

module load bioinfo/fastqc

echo "#### QC started `date` \n";

#fastqc --threads 20 --memory 10000 --outdir $outpath $inpath/Unknown_BK846-002L0001_good_1.fq.gz $inpath/Unknown_BK846-002L0002_good_2.fq.gz $inpath/Unknown_BK846-002L0004_good_1.fq.gz $inpath/Unknown_BK846-002L0005_good_2.fq.gz $inpath/Unknown_BK846-002L0001_good_2.fq.gz $inpath/Unknown_BK846-002L0003_good_1.fq.gz $inpath/Unknown_BK846-002L0004_good_2.fq.gz $inpath/Unknown_BK846-002L0006_good_1.fq.gz $inpath/Unknown_BK846-002L0002_good_1.fq.gz $inpath/Unknown_BK846-002L0003_good_2.fq.gz $inpath/Unknown_BK846-002L0005_good_1.fq.gz $inpath/Unknown_BK846-002L0006_good_2.fq.gz ;

fastqc --threads 20 --memory 10000 --outdir $outpath $inpath/Unknown_BK846-01L0013_good_1.fq.gz $inpath/Unknown_BK846-01L0014_good_2.fq.gz $inpath/Unknown_BK846-01L0016_good_1.fq.gz $inpath/Unknown_BK846-01L0017_good_2.fq.gz $inpath/Unknown_BK846-01L0013_good_2.fq.gz $inpath/Unknown_BK846-01L0015_good_1.fq.gz $inpath/Unknown_BK846-01L0016_good_2.fq.gz $inpath/Unknown_BK846-01L0018_good_1.fq.gz $inpath/Unknown_BK846-01L0014_good_1.fq.gz $inpath/Unknown_BK846-01L0015_good_2.fq.gz $inpath/Unknown_BK846-01L0017_good_1.fq.gz $inpath/Unknown_BK846-01L0018_good_2.fq.gz ;

echo "#### QC completed `date` \n";
wait



