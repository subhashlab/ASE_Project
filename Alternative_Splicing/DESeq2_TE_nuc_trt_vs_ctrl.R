library(DESeq2)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(readxl)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(tibble)
library(ggplot2)
library(ggrepel)

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


# Select nuclear samples with Control or HNRNPKsi
nuc_samples <- sample_info %>%
   filter(extraction == "nuclear", condition %in% c("Control", "HNRNPKsi"))
 
# Check
nuc_samples

nrow(nuc_samples)
#[1] 6
counts_nuc <- counts[, rownames(nuc_samples)]
dim(counts_nuc)   # sanity check
#[1] 4693511       6

counts_nuc_round <- round(counts_nuc)
# Had to use round function because TE counts are fractional.
dds_nuc <- DESeqDataSetFromMatrix(
  countData = counts_nuc_round,
  colData   = nuc_samples,
  design    = ~ condition
)
#converting counts to integer mode
#factor levels were dropped which had no samples


dds_nuc <- DESeq(dds_nuc)
res_nuc <- results(dds_nuc, contrast = c("condition", "HNRNPKsi", "Control"))
summary(res_nuc)
#out of 1651375 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)       : 5061, 0.31%
#LFC < 0 (down)     : 1745, 0.11%
#outliers [1]       : 0, 0%
#low counts [2]     : 1387787, 84%
#(mean count < 10)
#[1] see 'cooksCutoff' argument of ?results
#[2] see 'independentFiltering' argument of ?results

transcript_info <- read.csv("/iitjfs/home/2024rbb1003/projects/ASE_NB_project/RIPSeq/TE_annotation.csv", header = TRUE, stringsAsFactors = FALSE)

res_df <- as.data.frame(res_nuc) %>%                               
   rownames_to_column("transcript_id") %>%
   left_join(transcript_info, by = "transcript_id")

write.csv(res_df, "Nuc_HNRNPKsi_vs_Control_TE_results.csv", row.names = FALSE)
 

vsd <- vst(dds_nuc, blind = FALSE)
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


res_df <- as.data.frame(res_nuc) %>%                                                                                 
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
write.csv(res_sig, "nuc_HNRNPKsi_vs_Control_significant_0.05_1.csv", row.names = FALSE)

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
 
ggsave("Volcano_nuc_HNRNPKsi_vs_Control.pdf", volcano, width = 6, height = 6)

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

col_ha <- HeatmapAnnotation(condition = colData(dds_nuc)$condition,
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
 
pdf("Heatmap_nuc_HNRNPKsi_vs_Control_Top40_DEGs.pdf", width = 10, height = 12)
draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()
 

