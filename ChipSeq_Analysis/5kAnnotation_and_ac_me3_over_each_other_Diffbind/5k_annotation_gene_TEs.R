library(GenomicFeatures)
library(GenomicRanges)
library(rtracklayer)

# -----------------------------
# LOAD GTF
# -----------------------------

txdb_gene <- makeTxDbFromGFF(
  "/iitjfs/home/2025rbb1002/hg38/Homo_sapiens.GRCh38.113.chr.gtf",
  format = "gtf"
)

gtf <- import(
  "/iitjfs/home/2025rbb1002/hg38/Homo_sapiens.GRCh38.113.chr.gtf"
)

tx_map <- unique(
  as.data.frame(
    mcols(
      gtf[, c(
        "transcript_id",
        "gene_id",
        "gene_name",
        "gene_biotype"
      )]
    )
  )
)

# -----------------------------
# INPUT / OUTPUT FILES
# -----------------------------

peak_files <- c(
  "/iitjfs/home/2025rbb1002/M_ChipSeq/M_OP/batch2_t2p/DiffBind/summits_1000/K27ac/original_coordinates/oc/promoters/K27ac_all_down_peaks_merged_peaks.bed",
  "/iitjfs/home/2025rbb1002/M_ChipSeq/M_OP/batch2_t2p/DiffBind/summits_1000/K27ac/original_coordinates/oc/promoters/K27ac_all_up_peaks_merged_peaks.bed"
)

out_files <- c(
  "/iitjfs/home/2025rbb1002/M_ChipSeq/M_OP/batch2_t2p/DiffBind/summits_1000/K27ac/original_coordinates/oc/promoters/Annotated_diff_peaks_K27ac_sig_gene_down_peak5kb_TSS.csv",
  "/iitjfs/home/2025rbb1002/M_ChipSeq/M_OP/batch2_t2p/DiffBind/summits_1000/K27ac/original_coordinates/oc/promoters/Annotated_diff_peaks_K27ac_sig_gene_up_peak5kb_TSS.csv"
)

# -----------------------------
# TRANSCRIPTS
# -----------------------------

tx_gr <- transcripts(txdb_gene)

# Exact transcript TSSs
tss_exact <- promoters(
  tx_gr,
  upstream = 0,
  downstream = 1
)

# -----------------------------
# LOOP
# -----------------------------

for(i in seq_along(peak_files)) {
  
  bed <- read.table(
    peak_files[i],
    sep = "\t",
    header = FALSE,
    stringsAsFactors = FALSE
  )
  
  colnames(bed) <- c(
    "seqnames",
    "start",
    "end"
  )
  
  diff_peaks <- GRanges(
    seqnames = bed$seqnames,
    ranges = IRanges(
      start = bed$start,
      end = bed$end
    )
  )
  
  # Peak ±5 kb
  peak_5kb <- resize(
    diff_peaks,
    width = width(diff_peaks) + 10000,
    fix = "center"
  )
  
  # TSS overlap Peak ±5 kb
  hits <- findOverlaps(
    tss_exact,
    peak_5kb
  )
  
  tx_idx <- queryHits(hits)
  peak_idx <- subjectHits(hits)
  
  tss_pos <- start(tss_exact)[tx_idx]
  
  peak_start <- start(diff_peaks)[peak_idx]
  peak_end   <- end(diff_peaks)[peak_idx]
  
  distanceToTSS <- ifelse(
    tss_pos < peak_start,
    tss_pos - peak_start,
    ifelse(
      tss_pos > peak_end,
      tss_pos - peak_end,
      0
    )
  )
  
  final_anno <- data.frame(
    peak_index = peak_idx,
    seqnames = as.character(
      seqnames(diff_peaks)[peak_idx]
    ),
    start = peak_start,
    end = peak_end,
    transcript_id = mcols(tx_gr)$tx_name[tx_idx],
    tx_id = mcols(tx_gr)$tx_id[tx_idx],
    strand = as.character(
      strand(tx_gr)[tx_idx]
    ),
    distanceToTSS = distanceToTSS
  )
  
  final_anno <- merge(
    final_anno,
    tx_map,
    by = "transcript_id",
    all.x = TRUE
  )
  
  write.csv(
    final_anno,
    out_files[i],
    row.names = FALSE
  )
  
  cat(
    basename(out_files[i]),
    "\nOverlapping transcripts:", nrow(final_anno),
    "\nUnique transcripts:", length(unique(final_anno$transcript_id)),
    "\nUnique genes:", length(unique(final_anno$gene_id)),
    "\n\n"
  )
}
################################################################################
#TE Annotation 5k
################################################################################
library(GenomicRanges)
library(rtracklayer)

# -----------------------------
# LOAD TE ANNOTATION
# -----------------------------

te_gr <- import(
  "/iitjfs/home/2025rbb1002/hg38/GRCh38_Ensembl_rmsk_TE.gtf"
)

# -----------------------------
# INPUT / OUTPUT FILES
# -----------------------------

peak_files <- c(
  "/iitjfs/home/2025rbb1002/M_ChipSeq/M_OP/batch2_t2p/DiffBind/summits_1000/K27me3/original_coordinates/oc/promoters/K27me3_all_down_peaks_merged_peaks.bed",
  "/iitjfs/home/2025rbb1002/M_ChipSeq/M_OP/batch2_t2p/DiffBind/summits_1000/K27me3/original_coordinates/oc/promoters/K27me3_all_up_peaks_merged_peaks.bed"
)

out_files <- c(
  "/iitjfs/home/2025rbb1002/M_ChipSeq/M_OP/batch2_t2p/DiffBind/summits_1000/K27me3/original_coordinates/oc/promoters/Annotated_diff_peaks_K27me3_sig_TE_down_peak5kb_TSS.csv",
  "/iitjfs/home/2025rbb1002/M_ChipSeq/M_OP/batch2_t2p/DiffBind/summits_1000/K27me3/original_coordinates/oc/promoters/Annotated_diff_peaks_K27me3_sig_TE_up_peak5kb_TSS.csv"
)

# -----------------------------
# TE 5' END (strand-aware TSS)
# -----------------------------

te_tss <- resize(
  te_gr,
  width = 1,
  fix = "start"
)

# -----------------------------
# LOOP
# -----------------------------

for(i in seq_along(peak_files)) {
  
  bed <- read.table(
    peak_files[i],
    sep = "\t",
    header = FALSE,
    stringsAsFactors = FALSE
  )
  
  colnames(bed) <- c(
    "seqnames",
    "start",
    "end"
  )
  
  diff_peaks <- GRanges(
    seqnames = bed$seqnames,
    ranges = IRanges(
      start = bed$start,
      end = bed$end
    )
  )
  
  # -----------------------------
  # Peak ±5 kb
  # -----------------------------
  
  peak_5kb <- resize(
    diff_peaks,
    width = width(diff_peaks) + 10000,
    fix = "center"
  )
  
  # -----------------------------
  # TE TSS overlap Peak ±5 kb
  # -----------------------------
  
  hits_te <- findOverlaps(
    te_tss,
    peak_5kb
  )
  
  te_idx <- queryHits(hits_te)
  peak_idx <- subjectHits(hits_te)
  
  # -----------------------------
  # TE TSS position
  # -----------------------------
  
  te_tss_pos <- start(te_tss)[te_idx]
  
  peak_start <- start(diff_peaks)[peak_idx]
  peak_end   <- end(diff_peaks)[peak_idx]
  
  # -----------------------------
  # Distance from peak boundaries
  # -----------------------------
  
  distanceToTEstart <- ifelse(
    te_tss_pos < peak_start,
    te_tss_pos - peak_start,
    ifelse(
      te_tss_pos > peak_end,
      te_tss_pos - peak_end,
      0
    )
  )
  
  # -----------------------------
  # Final table
  # -----------------------------
  
  final_te_anno <- data.frame(
    peak_index = peak_idx,
    
    chr = as.character(
      seqnames(diff_peaks)[peak_idx]
    ),
    
    peak_start = peak_start,
    peak_end = peak_end,
    
    TE_chr = as.character(
      seqnames(te_gr)[te_idx]
    ),
    
    TE_start = start(te_gr)[te_idx],
    TE_end = end(te_gr)[te_idx],
    
    TE_strand = as.character(
      strand(te_gr)[te_idx]
    ),
    
    TE_id = mcols(te_gr)$gene_id[te_idx],
    
    transcript_id = mcols(te_gr)$transcript_id[te_idx],
    
    class_id = mcols(te_gr)$class_id[te_idx],
    
    distanceToTEstart = distanceToTEstart
  )
  
  write.csv(
    final_te_anno,
    out_files[i],
    row.names = FALSE
  )
  
  cat(
    basename(out_files[i]),
    "\nOverlaps:", nrow(final_te_anno),
    "\nUnique TEs:", length(unique(final_te_anno$TE_id)),
    "\n\n"
  )
}