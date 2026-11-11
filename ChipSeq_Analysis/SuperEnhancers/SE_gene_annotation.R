library(GenomicRanges)
library(rtracklayer)

# -----------------------------
# LOAD GTF
# -----------------------------
gtf <- import("/iitjfs/home/2025rbb1002/hg38/Homo_sapiens.GRCh38.113.chr.gtf")

genes <- gtf[gtf$type == "gene"]

# -----------------------------
# CREATE 1-bp TSS REGIONS
# -----------------------------
genes_tss <- promoters(genes, upstream = 0, downstream = 1)

# -----------------------------
# LOAD SE FILE
# BED must contain:
# chr start end SE_name
# -----------------------------
se_df <- read.table(
  "/iitjfs/home/2025rbb1002/M_ChipSeq/M_OP/batch2_t2p/DiffBind/summits_1000/K27ac/original_coordinates/oc/len/SE_MChip/diffbind_SE/Union/SE_downregulated.bed",
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
# (SE_start-500kb to SE_end+500kb)
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
# FIND GENES WHOSE TSS FALLS
# INSIDE THE ±500 kb WINDOW
# -----------------------------
hits <- findOverlaps(se_500kb, genes_tss)

# -----------------------------
# CREATE OUTPUT TABLE
# ONE ROW PER SE-GENE PAIR
# -----------------------------
tss_annotation <- data.frame(
  
  SE_name = mcols(se)$SE_ID[queryHits(hits)],
  
  SE_chr = as.character(seqnames(se))[queryHits(hits)],
  SE_start = start(se)[queryHits(hits)],
  SE_end = end(se)[queryHits(hits)],
  SE_width = width(se)[queryHits(hits)],
  
  window_start = start(se_500kb)[queryHits(hits)],
  window_end = end(se_500kb)[queryHits(hits)],
  
  gene_id = mcols(genes)$gene_id[subjectHits(hits)],
  gene_name = mcols(genes)$gene_name[subjectHits(hits)],
  gene_biotype = mcols(genes)$gene_biotype[subjectHits(hits)],
  
  gene_chr = as.character(seqnames(genes))[subjectHits(hits)],
  gene_start = start(genes)[subjectHits(hits)],
  gene_end = end(genes)[subjectHits(hits)],
  gene_width = width(genes)[subjectHits(hits)],
  
  strand = as.character(strand(genes))[subjectHits(hits)],
  
  TSS = start(genes_tss)[subjectHits(hits)]
)

# Remove duplicate SE-gene pairs if any
tss_annotation <- unique(tss_annotation)

# -----------------------------
# SAVE FULL ANNOTATION
# -----------------------------
write.csv(
  tss_annotation,
  "/iitjfs/home/2025rbb1002/M_ChipSeq/M_OP/batch2_t2p/DiffBind/summits_1000/K27ac/original_coordinates/oc/len/SE_MChip/diffbind_SE/Union/SE_down_TSS_within_500kb_annotation.csv",
  row.names = FALSE
)

# -----------------------------
# COUNT GENES PER SE
# -----------------------------
gene_counts <- aggregate(
  gene_id ~ SE_name,
  data = tss_annotation,
  FUN = function(x) length(unique(x))
)

colnames(gene_counts)[2] <- "genes_with_TSS_within_500kb"

write.csv(
  gene_counts,
  "/iitjfs/home/2025rbb1002/M_ChipSeq/M_OP/batch2_t2p/DiffBind/summits_1000/K27ac/original_coordinates/oc/len/SE_MChip/diffbind_SE/Union/SE_down_TSS_within_500kb_counts.csv",
  row.names = FALSE
)

# -----------------------------
# SUMMARY
# -----------------------------
cat("Number of SEs:", length(unique(tss_annotation$SE_name)), "\n")
cat("Number of unique genes:", length(unique(tss_annotation$gene_id)), "\n")
cat("Number of SE-gene pairs:", nrow(tss_annotation), "\n")
