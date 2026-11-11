python /iitjfs/home/2025rbb1002/rose/ROSE_main.py \
-g hg38 \
-i /iitjfs/home/2025rbb1002/M_ChipSeq/batch2_raw/Peaks/H3K27ac/control.gff \
-r /iitjfs/home/2025rbb1002/M_ChipSeq/batch2_raw/Unknown_CR300-001R0001.clean.bam /iitjfs/home/2025rbb1002/M_ChipSeq/batch2_raw/Unknown_CR300-001R0002.clean.bam \
-c /iitjfs/home/2025rbb1002/M_ChipSeq/batch2_raw/Unknown_CR300-001R0013.clean.bam /iitjfs/home/2025rbb1002/M_ChipSeq/batch2_raw/Unknown_CR300-001R0014.clean.bam \
-o rose_control

python /iitjfs/home/2025rbb1002/rose/ROSE_main.py \
-g hg38 \
-i /iitjfs/home/2025rbb1002/M_ChipSeq/batch2_raw/Peaks/H3K27ac/treat.gff \
-r /iitjfs/home/2025rbb1002/M_ChipSeq/batch2_raw/Unknown_CR300-001R0003.clean.bam /iitjfs/home/2025rbb1002/M_ChipSeq/batch2_raw/Unknown_CR300-001R0004.clean.bam /iitjfs/home/2025rbb1002/M_ChipSeq/batch2_raw/Unknown_CR300-001R0005.clean.bam /iitjfs/home/2025rbb1002/M_ChipSeq/batch2_raw/Unknown_CR300-001R0006.clean.bam \
-c /iitjfs/home/2025rbb1002/M_ChipSeq/batch2_raw/Unknown_CR300-001R0015.clean.bam /iitjfs/home/2025rbb1002/M_ChipSeq/batch2_raw/Unknown_CR300-001R0016.clean.bam /iitjfs/home/2025rbb1002/M_ChipSeq/batch2_raw/Unknown_CR300-001R0017.clean.bam /iitjfs/home/2025rbb1002/M_ChipSeq/batch2_raw/Unknown_CR300-001R0018.clean.bam \
-o rose_treat

cat /iitjfs/home/2025rbb1002/M_ChipSeq/M_OP/batch2_t2p/DiffBind/summits_1000/K27ac/original_coordinates/oc/len/SE_MChip/rose_control/control_SE.bed /iitjfs/home/2025rbb1002/M_ChipSeq/M_OP/batch2_t2p/DiffBind/summits_1000/K27ac/original_coordinates/oc/len/SE_MChip/rose_treat/treatment_SE.bed > /iitjfs/home/2025rbb1002/M_ChipSeq/M_OP/batch2_t2p/DiffBind/summits_1000/K27ac/original_coordinates/oc/len/SE_MChip/Combined_SE.bed

sort -k1,1 -k2,2n /iitjfs/home/2025rbb1002/M_ChipSeq/M_OP/batch2_t2p/DiffBind/summits_1000/K27ac/original_coordinates/oc/len/SE_MChip/Combined_SE.bed > /iitjfs/home/2025rbb1002/M_ChipSeq/M_OP/batch2_t2p/DiffBind/summits_1000/K27ac/original_coordinates/oc/len/SE_MChip/Combined_SE_sorted.bed

bedtools merge -i /iitjfs/home/2025rbb1002/M_ChipSeq/M_OP/batch2_t2p/DiffBind/summits_1000/K27ac/original_coordinates/oc/len/SE_MChip/Combined_SE_sorted.bed > /iitjfs/home/2025rbb1002/M_ChipSeq/M_OP/batch2_t2p/DiffBind/summits_1000/K27ac/original_coordinates/oc/len/SE_MChip/SE_union.bed

awk 'BEGIN{OFS="\t"} {print $1,$2,$3,"SE_"NR}' \
/iitjfs/home/2025rbb1002/M_ChipSeq/M_OP/batch2_t2p/DiffBind/summits_1000/K27ac/original_coordinates/oc/len/SE_MChip/diffbind_SE/Union/SE_union.bed \
> /iitjfs/home/2025rbb1002/M_ChipSeq/M_OP/batch2_t2p/DiffBind/summits_1000/K27ac/original_coordinates/oc/len/SE_MChip/diffbind_SE/Union/SE_union_named.bed
