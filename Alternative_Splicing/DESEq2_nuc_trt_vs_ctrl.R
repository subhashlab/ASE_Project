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

# Read counts file (keep rownames as gene IDs)
counts <- read.csv("/home/2024rbb1003/projects/ASE_NB_project/RNASeq/Quantification_gene/counts.csv", header = TRUE, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)

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
#[1] 78724     6
dds_nuc <- DESeqDataSetFromMatrix(
   countData = counts_nuc,
   colData   = nuc_samples,
   design    = ~ condition
 )
dds_nuc <- DESeq(dds_nuc)
res_nuc <- results(dds_nuc, contrast = c("condition", "HNRNPKsi", "Control"))
summary(res_nuc)

res_df <- as.data.frame(res_nuc) %>%
  rownames_to_column("gene_id")
 
write.csv(res_df, "nuc_HNRNPKsi_vs_Control_DESeq2.csv", row.names = FALSE)
 
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


res <- results(dds_nuc, contrast = c("condition", "HNRNPKsi", "Control"))
 
gene_info <- read_tsv("/iitjfs/home/2024rbb1003/projects/projectRseq/EndoOrganoids/scripts_combined/gene_info_hg38_113.tsv", show_col_types = FALSE)
res_df <- as.data.frame(res) %>%                                                                                 
   rownames_to_column("gene_id") %>%
   left_join(gene_info, by = "gene_id") %>%
   filter(!is.na(padj))
 
write.csv(res_df, "nuc_HNRNPKsi_vs_Control_results.csv", row.names = FALSE)

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
top_genes <- bind_rows(top_up, top_down)
 
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
     data = top_genes,
     aes(label = gene_name),
     size = 3,
     max.overlaps = 100
   ) +
   theme_bw() +
   labs(title = "HNRNPKsi vs Control", x = "Log2 Fold Change", y = "-Log10(Adjusted P-value)") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"), legend.title = element_blank())
 
ggsave("Volcano_nuc_HNRNPKsi_vs_Control.pdf", volcano, width = 6, height = 6)

vsd_mat <- assay(vsd)
sig_genes <- res_df %>%
  filter(padj < 0.05 & abs(log2FoldChange) >= 1)
 
deg_ids <- sig_genes$gene_id
mat_deg <- vsd_mat[deg_ids, , drop = FALSE]
rownames(mat_deg) <- sig_genes$gene_name
 
# Z-score scaling
mat_deg <- t(apply(mat_deg, 1, function(y) (y - mean(y)) / sd(y)^as.logical(sd(y))))
 
# Pick top 40 variable DEGs
gene_vars <- apply(mat_deg, 1, var)
top40_genes <- names(sort(gene_vars, decreasing = TRUE))[1:40]

col_ha <- HeatmapAnnotation(condition = colData(dds_nuc)$condition,
                              col = list(condition = c("Control"="#4472C4","HNRNPKsi"="#ed7d31")))
 
row_ha <- rowAnnotation(
   gene = anno_mark(
     at = which(rownames(mat_deg) %in% top40_genes),
     labels = rownames(mat_deg)[rownames(mat_deg) %in% top40_genes],
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
 
library(ggplot2)
 
# --- Set gene of interest ---
gene_of_interest <- "HDAC7"
gene_of_interest <- "EGR2"
gene_of_interest <- "FGL1"
gene_of_interest <- "PML"
# --- Map gene_name to gene_id ---
gene_id <- gene_info$gene_id[gene_info$gene_name == gene_of_interest]
if (length(gene_id) == 0) {
  stop(paste("Gene", gene_of_interest, "not found in gene_info mapping."))
}

# --- Match to rownames in vsd ---
gene_row <- gene_id[gene_id %in% rownames(vsd)]
if (length(gene_row) == 0) {
  stop(paste("Gene", gene_of_interest, "not found in vsd rownames."))
}
 
# --- Extract expression values ---
expr_values <- assay(vsd)[gene_row, ]
 
# --- Build plotting dataframe ---
df_gene <- data.frame(
  expression = as.numeric(expr_values),
  sample = colnames(vsd),
  group  = colData(vsd)$condition  # or $group depending on your DE design
)
 
# --- Create boxplot ---
p <- ggplot(df_gene, aes(x = group, y = expression, fill = group, color = group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.4) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.7) +
  theme_classic() +
  labs(title = gene_of_interest,
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
ggsave(paste0("expression_", gene_of_interest, "_boxplot_nuc.pdf"), p, width = 4, height = 4, dpi = 300)


