#module load bioinfo/gcc-13.1.0
#module load bioinfo/R-4.5.1
#export LD_LIBRARY_PATH=/iitjfs/apps/compiler/gcc-13.1.0/lib64
#R --vanilla


.libPaths(c(
   "/iitjfs/home/2024rbb1003/R/x86_64-pc-linux-gnu-library/4.5",
   "/iitjfs/apps/codes/R-4.5.1/lib64/R/library"
 ))

system("ldd ~/R/x86_64-pc-linux-gnu-library/4.5/Rcpp/libs/Rcpp.so | grep libstdc++")
library(DESeq2)
library(dplyr)
library(tibble)
library(readr)
library(dplyr)
library(readxl)

counts_raw <- read.table("/iitjfs/home/2024rbb1003/projects/ASE_NB_project/RNASeq/Differentially_expressed_TE_repeats/quantification/hg38_counts_TE.txt", header = TRUE, sep = "\t", check.names = FALSE)
head(counts_raw)
                                                                                  
# Use transcriptid as rownames
counts_raw <- counts_raw %>%
   column_to_rownames("Geneid")
head(counts_raw)
counts <- counts_raw[, -(1:5)]
head(counts)
colnames(counts) <- colnames(counts) %>%
  gsub(".*/", "", .) %>%
  gsub("Aligned.sortedByCoord.out.bam", "", .)
colnames(counts)
head(counts)
 
sample_map <- c(
  "Csi_S2_EKRN240050189-1A_22GGLCLT4" = "Csi_S2_EKRN240050189-1A",
  "Csi_S1_EKRN240050188-1A_22GGLCLT4" = "Csi_S1_EKRN240050188-1A",
  "ELAVL1si1_S1_EKRN240050190-1A_22GGLCLT4"  = "ELAVL1si1_S1_EKRN240050190-1A",
  "ELAVL1si1_S2_EKRN240050191-1A_22GGLCLT4"  = "ELAVL1si1_S2_EKRN240050191-1A",
  "ELAVL1si2_S1_EKRN240050192-1A_22GGLCLT4"  = "ELAVL1si2_S1_EKRN240050192-1A",
  "ELAVL1si2_S2_EKRN240050193-1A_22GGLCLT4"  = "ELAVL1si2_S2_EKRN240050193-1A"
)
 
old_names <- colnames(counts)
new_names <- old_names
idx <- match(old_names, names(sample_map))
 
new_names[!is.na(idx)] <- sample_map[idx[!is.na(idx)]]
new_names
colnames(counts) <- new_names
colnames(counts)
           
meta <- read_excel("/iitjfs/home/2024rbb1003/projects/ASE_NB_project/ASE_proj_sample_information_03102025.xlsx")
colnames(meta)                                                              

head(meta)
head(colnames(counts))
head(meta$sample_names)

sample_info <- meta%>%
  filter(sample_names %in% colnames(counts))
sample_info <- sample_info %>%
  arrange(match(sample_names, colnames(counts)))
all(sample_info$sample_names == colnames(counts))
#[1] TRUE
colnames(counts)

sample_info$sample_names
head(sample_info)
# Save sample_info
write.csv(sample_info, "sample_info.csv", row.names = FALSE)
 
# Save counts
write.csv(counts, "counts_te.csv", row.names = TRUE)  # keep rownames as transcript IDs





library(DESeq2)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(readxl)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(tibble)
# Read the sample_info file
sample_info <- read.csv("/iitjfs/home/2024rbb1003/projects/ASE_NB_project/RNASeq/Differentially_expressed_TE_repeats/results/sample_info.csv", header = TRUE, stringsAsFactors = FALSE)

# Read counts file (keep rownames as transcript IDs)
counts <- read.csv("/iitjfs/home/2024rbb1003/projects/ASE_NB_project/RNASeq/Differentially_expressed_TE_repeats/results/counts_te.csv", header = TRUE, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)

sample_info <- as.data.frame(sample_info)
rownames(sample_info) <- sample_info$sample_names

setdiff(colnames(counts), rownames(sample_info))
#character(0)

sample_info <- sample_info %>%
  mutate(
    condition  = factor(condition),
    treatment  = factor(treatment),
    cellline   = factor(cellline),
    extraction = factor(extraction),
    batch      = factor(batch)
  )
head(sample_info)


# Select cytoplasmic samples with Control or HNRNPKsi
cyto_samples <- sample_info %>%
   filter(extraction == "cytoplasmic", condition %in% c("Control", "HNRNPKsi"))
 
# Check
cyto_samples

nrow(cyto_samples)
#[1] 6
counts_cyto <- counts[, rownames(cyto_samples)]
dim(counts_cyto)   # sanity check
#[1] 4693511       6

counts_cyto_round <- round(counts_cyto)
# Had to use round function because TE counts are fractional.
dds_cyto <- DESeqDataSetFromMatrix(
  countData = counts_cyto_round,
  colData   = cyto_samples,
  design    = ~ condition
)
#converting counts to integer mode
#factor levels were dropped which had no samples

dds_cyto <- DESeq(dds_cyto)
res_cyto <- results(dds_cyto, contrast = c("condition", "HNRNPKsi", "Control"))
summary(res_cyto)
#out of 544885 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)       : 2183, 0.4%
#LFC < 0 (down)     : 1382, 0.25%
#outliers [1]       : 0, 0%
#low counts [2]     : 308646, 57%
#(mean count < 1)
#[1] see 'cooksCutoff' argument of ?results
#[2] see 'independentFiltering' argument of ?results

transcript_info <- read.csv("/iitjfs/home/2024rbb1003/projects/ASE_NB_project/RIPSeq/TE_annotation.csv", header = TRUE, stringsAsFactors = FALSE)

res_df <- as.data.frame(res_cyto) %>%                               
   rownames_to_column("transcript_id") %>%
   left_join(transcript_info, by = "transcript_id")

write.csv(res_df, "Cyto_HNRNPKsi_vs_Control_TE_results.csv", row.names = FALSE)
 

vsd <- vst(dds_cyto, blind = FALSE)
pcaData <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)

percentVar <- round(100 * attr(pcaData, "percentVar"))
pcaData$sample_name <- rownames(pcaData)
 
# PCA WITH labels
p_with <- ggplot(pcaData, aes(PC1, PC2, color = condition)) +
   geom_point(size = 3) +
   geom_text_repel(aes(label = sample_name), size = 3, max.overlaps = Inf) +
   xlab(paste0("PC1: ", percentVar[1], "% variance")) +
   ylab(paste0("PC2: ", percentVar[2], "% variance")) +
   theme_bw(base_size = 12) +
   ggtitle("PCA (batch-corrected): HNRNPKsi vs Control")

# PCA WITHOUT labels
p_no <- ggplot(pcaData, aes(PC1, PC2, color = condition)) +
   geom_point(size = 3) +
   xlab(paste0("PC1: ", percentVar[1], "% variance")) +
   ylab(paste0("PC2: ", percentVar[2], "% variance")) +
   theme_bw(base_size = 12) +
   ggtitle("PCA (batch-corrected): HNRNPKsi vs Control")

# Save PCA plots
ggsave("PCA_HNRNPKsi_vs_Control_with_labels.pdf", p_with, width = 7, height = 6)
ggsave("PCA_HNRNPKsi_vs_Control_no_labels.pdf", p_no, width = 7, height = 6)

res_df <- as.data.frame(res_cyto) %>%                                                                                 
   rownames_to_column("transcript_id") %>%
   left_join(transcript_info, by = "transcript_id") %>%
   filter(!is.na(padj))

res_df <- res_df %>%
   mutate(category = case_when(
     padj < 0.05 & log2FoldChange >= 1 ~ "Upregulated",
     padj < 0.05 & log2FoldChange <= -1 ~ "Downregulated",
     TRUE ~ "Not Significant"
  ))
 
res_sig <- res_df %>% filter(category %in% c("Upregulated", "Downregulated"))

# Save to file (CSV or TSV)
write.csv(res_sig, "Cyto_HNRNPKsi_vs_Control_significant_0.05_1.csv", row.names = FALSE)

top_up <- res_df %>% filter(category == "Upregulated") %>% arrange(desc(log2FoldChange)) %>% head(10)
top_down <- res_df %>% filter(category == "Downregulated") %>% arrange(log2FoldChange) %>% head(10)
top_transcripts <- bind_rows(top_up, top_down)
 
volcano <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj))) +
   geom_point(aes(color = category), alpha = 0.6, size = 2) +
   scale_color_manual(values = c(
     "Upregulated" = "#ed7d31",
     "Downregulated" = "#0b5394",
     "Not Significant" = "gray70"
   )) +
   geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "gray50") +
   geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray50") +
   geom_text_repel(
     data = top_transcripts,
     aes(label = transcript_id),
     size = 3,
     max.overlaps = 100
   ) +
   theme_bw() +
   labs(title = "HNRNPKsi vs Control", x = "Log2 Fold Change", y = "-Log10(Adjusted P-value)") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"), legend.title = element_blank())
 
ggsave("Volcano_Cyto_HNRNPKsi_vs_Control_2.pdf", volcano, width = 6, height = 6)


vsd_mat <- assay(vsd)
sig_transcripts <- res_df %>%
  filter(padj < 0.05 & abs(log2FoldChange) >= 1)
 
deg_ids <- sig_transcripts$transcript_id
mat_deg <- vsd_mat[deg_ids, , drop = FALSE]

# Z-score scaling
mat_deg <- t(apply(mat_deg, 1, function(y) {
  if (sd(y) == 0) return(rep(0, length(y)))
  (y - mean(y)) / sd(y)
}))

 
# Pick top 40 variable DEGs
transcript_vars <- apply(mat_deg, 1, var)
top40_transcripts <- names(sort(transcript_vars, decreasing = TRUE))[1:40]

col_ha <- HeatmapAnnotation(condition = colData(dds_cyto)$condition,
                              col = list(condition = c("Control"="#4472C4","HNRNPKsi"="#ed7d31")))
 
row_ha <- rowAnnotation(
  transcript = anno_mark(
    at = which(rownames(mat_deg) %in% top40_transcripts),
    labels = top40_transcripts,
    labels_gp = gpar(fontsize = 8, fontface = "italic")
  )
)

 
col_fun <- colorRamp2(c(min(mat_deg), 0, max(mat_deg)), c("#2166AC","#F0F0F0","#E66101"))
 
ht <- Heatmap(mat_deg,
               name = "Z-Score",
               col = col_fun,
               top_annotation = col_ha,
               right_annotation = row_ha,
               show_row_names = FALSE,
               cluster_rows = TRUE,
               cluster_columns = TRUE,
               column_names_gp = gpar(fontsize = 6),
               row_names_gp = gpar(fontsize = 6),
               column_title = "Heatmap: Top DEGs HNRNPKsi vs Control")
 
pdf("Heatmap_Cyto_HNRNPKsi_vs_Control_Top40_DEGs.pdf", width = 10, height = 12)
draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()
 

