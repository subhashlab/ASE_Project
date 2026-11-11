library(GenomicRanges)
library(rtracklayer)

# -----------------------------
# LOAD TE GTF
# -----------------------------
te_gtf <- import("/iitjfs/home/2025rbb1002/hg38/GRCh38_Ensembl_rmsk_TE.gtf")

# Keep TE entries
tes <- te_gtf[te_gtf$type == "exon"]

# -----------------------------
# LOAD SE FILE
# BED: chr start end SE_name
# -----------------------------
se_df <- read.table(
  "/iitjfs/home/2025rbb1002/M_ChipSeq/M_OP/batch2_t2p/DiffBind/summits_1000/K27ac/original_coordinates/oc/len/SE_MChip/diffbind_SE/Union/SE_upregulated.bed",
  sep = ",",
  header = TRUE,
  stringsAsFactors = FALSE
)

#colnames(se_df) <- c("chr","start","end","SE_name")

se <- GRanges(
  seqnames = se_df$chr,
  ranges = IRanges(
    start = se_df$start + 1,  # BED is 0-based
    end = se_df$end
  )
)

mcols(se)$SE_ID <- se_df$SE_name


# -----------------------------
# CREATE SE ±500 kb WINDOW
# -----------------------------
se_500kb <- GRanges(
  seqnames = seqnames(se),
  ranges = IRanges(
    start = pmax(1, start(se) - 500000),
    end = end(se) + 500000
  )
)

mcols(se_500kb)$SE_ID <- mcols(se)$SE_ID

# -----------------------------
# FIND TEs OVERLAPPING WINDOW
# -----------------------------
hits <- findOverlaps(se_500kb, tes)

# -----------------------------
# CREATE OUTPUT TABLE
# ONE ROW PER SE–TE PAIR
# -----------------------------
te_annotation <- data.frame(

  SE_name = mcols(se)$SE_ID[queryHits(hits)],

  SE_chr = as.character(seqnames(se))[queryHits(hits)],
  SE_start = start(se)[queryHits(hits)],
  SE_end = end(se)[queryHits(hits)],
  SE_width = width(se)[queryHits(hits)],

  window_start = start(se_500kb)[queryHits(hits)],
  window_end = end(se_500kb)[queryHits(hits)],

  transcript_id = mcols(tes)$transcript_id[subjectHits(hits)],
  gene_id = mcols(tes)$gene_id[subjectHits(hits)],
  family_id = mcols(tes)$family_id[subjectHits(hits)],
  class_id = mcols(tes)$class_id[subjectHits(hits)],

  TE_chr = as.character(seqnames(tes))[subjectHits(hits)],
  TE_start = start(tes)[subjectHits(hits)],
  TE_end = end(tes)[subjectHits(hits)],
  TE_width = width(tes)[subjectHits(hits)],
  strand = as.character(strand(tes))[subjectHits(hits)]
)

# Remove duplicate SE-TE pairs
te_annotation <- unique(te_annotation)

# -----------------------------
# SAVE FULL ANNOTATION
# -----------------------------
write.csv(
  te_annotation,
  "/iitjfs/home/2025rbb1002/M_ChipSeq/M_OP/batch2_t2p/DiffBind/summits_1000/K27ac/original_coordinates/oc/len/SE_MChip/diffbind_SE/Union/SE_up_TE_within_500kb_annotation.csv",
  row.names = FALSE
)

# -----------------------------
# COUNT TEs PER SE
# -----------------------------
te_counts <- aggregate(
  transcript_id ~ SE_name,
  data = te_annotation,
  FUN = function(x) length(unique(x))
)

colnames(te_counts)[2] <- "TEs_within_500kb"

write.csv(
  te_counts,
  "/iitjfs/home/2025rbb1002/M_ChipSeq/M_OP/batch2_t2p/DiffBind/summits_1000/K27ac/original_coordinates/oc/len/SE_MChip/diffbind_SE/Union/SE_up_TE_within_500kb_counts.csv",
  row.names = FALSE
)

# -----------------------------
# FAMILY COUNTS
# -----------------------------
family_counts <- aggregate(
  family_id ~ SE_name,
  data = te_annotation,
  FUN = function(x) length(unique(x))
)

colnames(family_counts)[2] <- "TE_families_within_500kb"

write.csv(
  family_counts,
  "/iitjfs/home/2025rbb1002/M_ChipSeq/M_OP/batch2_t2p/DiffBind/summits_1000/K27ac/original_coordinates/oc/len/SE_MChip/diffbind_SE/Union/SE_up_TE_family_counts.csv",
  row.names = FALSE
)

# -----------------------------
# CLASS COUNTS
# -----------------------------
class_counts <- aggregate(
  class_id ~ SE_name,
  data = te_annotation,
  FUN = function(x) length(unique(x))
)

colnames(class_counts)[2] <- "TE_classes_within_500kb"

write.csv(
  class_counts,
  "/iitjfs/home/2025rbb1002/M_ChipSeq/M_OP/batch2_t2p/DiffBind/summits_1000/K27ac/original_coordinates/oc/len/SE_MChip/diffbind_SE/Union/SE_up_TE_class_counts.csv",
  row.names = FALSE
)

# -----------------------------
# SUMMARY
# -----------------------------
cat("Number of SEs:", length(unique(te_annotation$SE_name)), "\n")
cat("Number of unique TEs:", length(unique(te_annotation$transcript_id)), "\n")
cat("Number of unique families:", length(unique(te_annotation$family_id)), "\n")
cat("Number of unique classes:", length(unique(te_annotation$class_id)), "\n")
cat("Number of SE-TE pairs:", nrow(te_annotation), "\n")
