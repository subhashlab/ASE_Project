#HDAC
computeMatrix scale-regions \
    -S ac_avg_cntrl.bw \
    -R HDAC_tot_ri_2.bed \
    --missingDataAsZero \
    -o HDAC_tot_ri_ac_avg_cntrl_matrix.gz
    
computeMatrix scale-regions \
    -S ac_avg_hn.bw \
    -R HDAC_tot_ri_2.bed \
    --missingDataAsZero \
    -o HDAC_tot_ri_ac_avg_hn_matrix.gz
    
#Zoom-In
computeMatrix reference-point \
    --referencePoint center \
    -S /home/2025rbb1002/M_ChipSeq/M_OP/batch2_t2p/DiffBind/summits_1000/integration/ASE/ASE_heatmaps/plots/ac_avg_cntrl.bw \
    -R /home/2025rbb1002/M_ChipSeq/M_OP/batch2_t2p/DiffBind/summits_1000/integration/ASE/ASE_heatmaps/plots/tot_ri_ZoomIn/HDAC_tot.bed \
    -b 400 \
    -a 400 \
    --missingDataAsZero \
    -o /home/2025rbb1002/M_ChipSeq/M_OP/batch2_t2p/DiffBind/summits_1000/integration/ASE/ASE_heatmaps/plots/tot_ri_ZoomIn/HDAC_tot_ri_ac_avg_cntrl_matrix.gz
    
computeMatrix reference-point \
    --referencePoint center \
    -S /home/2025rbb1002/M_ChipSeq/M_OP/batch2_t2p/DiffBind/summits_1000/integration/ASE/ASE_heatmaps/plots/ac_avg_hn.bw \
    -R /home/2025rbb1002/M_ChipSeq/M_OP/batch2_t2p/DiffBind/summits_1000/integration/ASE/ASE_heatmaps/plots/tot_ri_ZoomIn/HDAC_tot.bed \
    -b 400 \
    -a 400 \
    --missingDataAsZero \
    -o /home/2025rbb1002/M_ChipSeq/M_OP/batch2_t2p/DiffBind/summits_1000/integration/ASE/ASE_heatmaps/plots/tot_ri_ZoomIn/HDAC_tot_ri_ac_avg_hn_matrix.gz
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
library(data.table)
library(ComplexHeatmap)
library(circlize)

# =========================================================
# FILES
# =========================================================

file1 <- "tot_ac_avg_cntrl_matrix_3.gz"
file2 <- "tot_ac_avg_hn_matrix_3.gz"

# =========================================================
# FUNCTION TO READ MATRIX
# =========================================================

read_matrix <- function(f) {

  x <- fread(
    cmd = paste("zcat", f),
    skip = 1,
    header = FALSE
  )

  # First 6 columns = genomic information
  info <- x[, 1:6]
  region_id <- x$V4
  # Columns 7 onward = signal
  signal <- as.matrix(x[, 7:ncol(x)])
  mode(signal) <- "numeric"

  # Row-wise Z-score
  signal <- t(scale(t(signal)))

  # Remove zero-variance / NA rows
  keep <- complete.cases(signal)

  info <- info[keep, ]
  signal <- signal[keep, ]

  # Clip
  signal[signal > 4] <- 4
  signal[signal < -4] <- -4

  return(
  list(
    info = info,
    signal = signal,
    id = region_id[keep]
  )
)
}


# =========================================================
# READ BOTH MATRICES
# =========================================================

m1 <- read_matrix(file1)
m2 <- read_matrix(file2)


# =========================================================
# CREATE UNIQUE REGION ID
# =========================================================



m1$id <- m1$id
m2$id <- m2$id

# =========================================================
# KEEP ONLY REGIONS PRESENT IN BOTH MATRICES
# =========================================================

common <- intersect(m1$id, m2$id)

cat("Matrix 1 regions:", length(m1$id), "\n")
cat("Matrix 2 regions:", length(m2$id), "\n")
cat("Common regions  :", length(common), "\n")


m1_keep <- match(common, m1$id)
m2_keep <- match(common, m2$id)

m1_signal <- m1$signal[m1_keep, ]
m2_signal <- m2$signal[m2_keep, ]

m1_info <- m1$info[m1_keep, ]


# =========================================================
# DETERMINE ORDER FROM MATRIX 1
# =========================================================

max_pos <- apply(
  m1_signal,
  1,
  which.max
)

gene_order <- order(max_pos)


# =========================================================
# APPLY SAME ORDER TO BOTH
# =========================================================

m1_signal <- m1_signal[gene_order, ]

m2_signal <- m2_signal[gene_order, ]

m1_info <- m1_info[gene_order, ]
rownames(m1_signal) <- m1$id[m1_keep][gene_order]
rownames(m2_signal) <- m2$id[m2_keep][gene_order]

# =========================================================
# CHECK
# =========================================================

head(m1_info)

head(max_pos[gene_order])


# =========================================================
# HEATMAP 1
# =========================================================

pdf(
  "tot_ac_avg_cntrl_matrix_3.pdf",
  width = 5,
  height = 8
)

draw(
  Heatmap(
    m1_signal,
    name = "Z-score",
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_names = TRUE,
    show_column_names = FALSE,
    use_raster = TRUE,
    col = colorRamp2(
      c(-4, 0, 4),
      c("blue", "white", "red")
    )
  )
)

dev.off()


# =========================================================
# HEATMAP 2
# =========================================================

pdf(
  "tot_ac_avg_hn_matrix_3.pdf",
  width = 5,
  height = 8
)

draw(
  Heatmap(
    m2_signal,
    name = "Z-score",
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_names = TRUE,
    show_column_names = FALSE,
    use_raster = TRUE,
    col = colorRamp2(
      c(-4, 0, 4),
      c("blue", "white", "red")
    )
  )
)

dev.off()
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
library(data.table)
library(ComplexHeatmap)
library(circlize)

# =========================================================
# FILES
# =========================================================

file1 <- "tot_ac_avg_cntrl_matrix_3.gz"
file2 <- "tot_ac_avg_hn_matrix_3.gz"


# =========================================================
# READ ORIGINAL SIGNAL MATRIX
# =========================================================

read_matrix <- function(f) {

  x <- fread(
    cmd = paste("zcat", f),
    skip = 1,
    header = FALSE
  )

  # First 6 columns = genomic information
  info <- x[, 1:6]
  region_id <- x$V4
  # Columns 7 onward = original signal
  signal <- as.matrix(x[, 7:ncol(x)])
  mode(signal) <- "numeric"

  list(
    info = info,
    signal = signal,
    id = region_id
  )
}


# =========================================================
# READ BOTH MATRICES
# =========================================================

m1 <- read_matrix(file1)
m2 <- read_matrix(file2)


# =========================================================
# CREATE REGION IDs
# =========================================================


m1$id <- m1$id
m2$id <- m2$id

# =========================================================
# MATCH COMMON REGIONS
# =========================================================

common <- intersect(m1$id, m2$id)

cat("Matrix 1 regions:", length(m1$id), "\n")
cat("Matrix 2 regions:", length(m2$id), "\n")
cat("Common regions  :", length(common), "\n")

m1_keep <- match(common, m1$id)
m2_keep <- match(common, m2$id)

m1_signal <- m1$signal[m1_keep, ]
m2_signal <- m2$signal[m2_keep, ]

m1_info <- m1$info[m1_keep, ]


# =========================================================
# REMOVE NA / ZERO-VARIANCE REGIONS
# =========================================================

sd1 <- apply(m1_signal, 1, sd)
sd2 <- apply(m2_signal, 1, sd)

keep <- (
  complete.cases(m1_signal) &
  complete.cases(m2_signal) &
  is.finite(sd1) &
  is.finite(sd2) &
  sd1 > 0 &
  sd2 > 0
)

cat("Regions after filtering:", sum(keep), "\n")

m1_signal <- m1_signal[keep, ]
m2_signal <- m2_signal[keep, ]
m1_info <- m1_info[keep, ]


# =========================================================
# DETERMINE ORDER FROM CONTROL
# =========================================================

max_pos <- apply(
  m1_signal,
  1,
  which.max
)

gene_order <- order(max_pos)


# =========================================================
# APPLY SAME ORDER TO BOTH
# =========================================================

m1_signal <- m1_signal[gene_order, ]
m2_signal <- m2_signal[gene_order, ]
m1_info <- m1_info[gene_order, ]
rownames(m1_signal) <- m1$id[m1_keep][keep][gene_order]
rownames(m2_signal) <- m2$id[m2_keep][keep][gene_order]

# =========================================================
# BIN-BY-BIN LOG2 FOLD CHANGE
# =========================================================

pseudocount <- 1

log2FC <- log2(
  (m2_signal + pseudocount) /
  (m1_signal + pseudocount)
)


# =========================================================
# CLIP LOG2FC
# =========================================================

log2FC[log2FC > 0.5] <- 0.5
log2FC[log2FC < -0.5] <- -0.5


# =========================================================
# PLOT
# =========================================================

pdf(
  "tot_ac_avg_HN_vs_Control_log2FC_3.pdf",
  width = 5,
  height = 8
)

draw(
  Heatmap(
    log2FC,
    name = "log2FC",
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_names = TRUE,
    show_column_names = FALSE,
    use_raster = TRUE,

    col = colorRamp2(
      c(-0.5, 0, 0.5),
      c("blue", "white", "red")
    )
  )
)

dev.off()
