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
sample_info <- read.csv("/home/2024rbb1003/projects/ASE_NB_project/RNASeq/Quantification_gene/sample_info.csv", header = TRUE, stringsAsFactors = FALSE)

# Read counts file (keep rownames as transcript IDs)
counts <- read.csv("/home/2024rbb1003/projects/ASE_NB_project/RNASeq/Quantification_gene/counts_transcript_level.csv", header = TRUE, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)

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


# Select total samples with Control or HNRNPKsi
Tot_samples <- sample_info %>%
   filter(extraction == "total", condition %in% c("Control", "HNRNPKsi"), batch %in% c("b1"))
 
# Check
Tot_samples

nrow(Tot_samples)
#[1] 6
counts_Tot <- counts[, rownames(Tot_samples)]
dim(counts_Tot)   # sanity check
#[1] 78724     6
dds_Tot <- DESeqDataSetFromMatrix(
   countData = counts_Tot,
   colData   = Tot_samples,
   design    = ~ condition
 )
dds_Tot <- DESeq(dds_Tot)
res_Tot <- results(dds_Tot, contrast = c("condition", "HNRNPKsi", "Control"))
summary(res_Tot)

transcript_info <- read_tsv("/iitjfs/home/2024rbb1003/projects/ASE_NB_project/transcript_gene_info.tsv", show_col_types = FALSE)

res_df <- as.data.frame(res_Tot) %>%                               
   rownames_to_column("transcript_id") %>%
   left_join(transcript_info, by = "transcript_id")

write.csv(res_df, "Tot_HNRNPKsi_vs_Control_results.csv", row.names = FALSE)
 

vsd <- vst(dds_Tot, blind = FALSE)
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


res <- results(dds_Tot, contrast = c("condition", "HNRNPKsi", "Control"))
 
transcript_info <- read_tsv("/iitjfs/home/2024rbb1003/projects/ASE_NB_project/transcript_gene_info.tsv", show_col_types = FALSE)
res_df <- as.data.frame(res) %>%                                                                                 
   rownames_to_column("transcript_id") %>%
   left_join(transcript_info, by = "transcript_id") %>%
   filter(!is.na(padj))
 
write.csv(res_df, "Tot_HNRNPKsi_vs_Control_results.csv", row.names = FALSE)

res_df <- res_df %>%
   mutate(category = case_when(
     padj < 0.05 & log2FoldChange >= 1 ~ "Upregulated",
     padj < 0.05 & log2FoldChange <= -1 ~ "Downregulated",
     TRUE ~ "Not Significant"
  ))
 
res_sig <- res_df %>% filter(category %in% c("Upregulated", "Downregulated"))

# Save to file (CSV or TSV)
write.csv(res_sig, "Tot_HNRNPKsi_vs_Control_significant_0.05_1.csv", row.names = FALSE)

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
     aes(label = transcript_name),
     size = 3,
     max.overlaps = 100
   ) +
   theme_bw() +
   labs(title = "HNRNPKsi vs Control", x = "Log2 Fold Change", y = "-Log10(Adjusted P-value)") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"), legend.title = element_blank())
 
ggsave("Volcano_Tot_HNRNPKsi_vs_Control.pdf", volcano, width = 6, height = 6)

vsd_mat <- assay(vsd)
sig_transcripts <- res_df %>%
  filter(padj < 0.05 & abs(log2FoldChange) >= 1)
 
deg_ids <- sig_transcripts$transcript_id
mat_deg <- vsd_mat[deg_ids, , drop = FALSE]
tx_lookup <- sig_transcripts %>%
  distinct(transcript_id, transcript_name) %>%
  column_to_rownames("transcript_id")

# Z-score scaling
mat_deg <- t(apply(mat_deg, 1, function(y) {
  if (sd(y) == 0) return(rep(0, length(y)))
  (y - mean(y)) / sd(y)
}))

 
# Pick top 40 variable DEGs
transcript_vars <- apply(mat_deg, 1, var)
top40_transcripts <- names(sort(transcript_vars, decreasing = TRUE))[1:40]

col_ha <- HeatmapAnnotation(condition = colData(dds_Tot)$condition,
                              col = list(condition = c("Control"="#4472C4","HNRNPKsi"="#ed7d31")))
 
row_ha <- rowAnnotation(
  transcript = anno_mark(
    at = which(rownames(mat_deg) %in% top40_transcripts),
    labels = tx_lookup[ top40_transcripts, "transcript_name" ],
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
 
pdf("Heatmap_Tot_HNRNPKsi_vs_Control_Top40_DEGs.pdf", width = 10, height = 12)
draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()
 
library(ggplot2)
 
# --- Set transcript of interest ---
transcript_of_interest <- "HDAC7-217"
transcript_of_interest <- "HDAC7-203"
transcript_of_interest <- "HDAC7-225"
transcript_of_interest <- "HDAC7-214"
transcript_of_interest <- "HDAC7-216"
transcript_of_interest <- "PML-210"
transcript_of_interest <- "PML-212"

# --- Map transcript_name to transcript_id ---
transcript_id <- transcript_info$transcript_id[transcript_info$transcript_name == transcript_of_interest]
if (length(transcript_id) == 0) {
  stop(paste("transcript", transcript_of_interest, "not found in transcript_info mapping."))
}

# --- Match to rownames in vsd ---
transcript_row <- transcript_id[transcript_id %in% rownames(vsd)]
if (length(transcript_row) == 0) {
  stop(paste("transcript", transcript_of_interest, "not found in vsd rownames."))
}
 
# --- Extract expression values ---
expr_values <- assay(vsd)[transcript_row, ]
 
# --- Build plotting dataframe ---
df_transcript <- data.frame(
  expression = as.numeric(expr_values),
  sample = colnames(vsd),
  group  = colData(vsd)$condition  # or $group depending on your DE design
)
 
# --- Create boxplot ---
p <- ggplot(df_transcript, aes(x = group, y = expression, fill = group, color = group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.4) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.7) +
  theme_classic() +
  labs(title = transcript_of_interest,
       x = "Group",
       y = "Normalized Expression (VST)") +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8)
  ) +
  scale_fill_manual(
   values = c("Control" = "#4472C4", "HNRNPKsi" = "#ed7d31")
  ) +
  scale_color_manual(
    values = c("Control" = "#4472C4", "HNRNPKsi" = "#ed7d31")
  ) +
  scale_x_discrete(
    labels = c("Control" = "Control", "HNRNPKsi" = "HNRNPKsi")
  )
 
# --- Save plot ---
ggsave(paste0("expression_", transcript_of_interest, "_boxplot_Tot.pdf"), p, width = 4, height = 4, dpi = 300)


