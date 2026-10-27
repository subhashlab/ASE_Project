
# ============================================================
# Load packages
# ============================================================

library(VennDiagram)
library(grid)

eclip1_file <- "/iitjfs/home/2024rbb1003/projects/ASE_NB_project/eCLIP/results/eclip_peaks1_annotated.csv"
eclip2_file <- "/iitjfs/home/2024rbb1003/projects/ASE_NB_project/eCLIP/results/eclip_peaks2_annotated.csv"


# ============================================================
# 2. Read eCLIP annotation files
# ============================================================

eclip1 <- read.csv(eclip1_file, header = TRUE, stringsAsFactors = FALSE)
eclip2 <- read.csv(eclip2_file, header = TRUE, stringsAsFactors = FALSE)

cat("eCLIP 1 genes:", nrow(eclip1), "\n")
cat("eCLIP 2 genes:", nrow(eclip2), "\n")


# ============================================================
# 3. Get common genes between the two eCLIP datasets
# ============================================================

eclip1_genes <- unique(na.omit(eclip1$SYMBOL))
eclip2_genes <- unique(na.omit(eclip2$SYMBOL))

common_eclip_genes <- intersect(
  eclip1_genes,
  eclip2_genes
)

cat("Unique genes in eCLIP 1:", length(eclip1_genes), "\n")
cat("Unique genes in eCLIP 2:", length(eclip2_genes), "\n")
cat("Common eCLIP genes:", length(common_eclip_genes), "\n")

parent_dir <- "/iitjfs/home/2024rbb1003/projects/ASE_NB_project/RNASeq/Alternative_splice_sites/rmats_results/rmats_filtered_results/"

folders <- c(
  "rmats_cyto_trt_vs_ctrl",
  "rmats_Nuc_trt_vs_ctrl",
  "rmats_Total_trt_vs_ctrl"
)


# ============================================================
# Read significant rMATS genes
# ============================================================

pattern <- "_JCEC_sig_FDR0.05_dPSI0.1.txt$"

rmats_gene_sets <- list()

for (folder in folders) {
  
  folder_path <- file.path(parent_dir, folder)
  
  files <- list.files(
    folder_path,
    pattern = pattern,
    full.names = TRUE
  )
  
  genes <- character(0)
  
  for (file in files) {
    
    dat <- read.delim(
      file,
      header = TRUE,
      sep = "\t",
      stringsAsFactors = FALSE
    )
    
    if (!"geneSymbol" %in% colnames(dat)) {
      warning(
        "geneSymbol column not found in: ",
        basename(file)
      )
      next
    }
    
    genes <- c(
      genes,
      dat$geneSymbol
    )
  }
  
  genes <- unique(
    na.omit(
      genes[genes != ""]
    )
  )
  
  rmats_gene_sets[[folder]] <- genes
  
  cat(
    folder,
    ":",
    length(genes),
    "unique rMATS genes\n"
  )
}


# ============================================================
# Assign individual gene sets
# ============================================================

cyto_genes <- rmats_gene_sets[["rmats_cyto_trt_vs_ctrl"]]

nuc_genes <- rmats_gene_sets[["rmats_Nuc_trt_vs_ctrl"]]

total_genes <- rmats_gene_sets[["rmats_Total_trt_vs_ctrl"]]


# ============================================================
# Print gene-set sizes
# ============================================================

cat("Cyto genes:", length(cyto_genes), "\n")
cat("Nuc genes:", length(nuc_genes), "\n")
cat("Total genes:", length(total_genes), "\n")
cat("Common eCLIP genes:", length(common_eclip_genes), "\n")


# ============================================================
# Calculate overlaps
# ============================================================

cat(
  "Cyto ∩ eCLIP:",
  length(intersect(cyto_genes, common_eclip_genes)),
  "\n"
)

cat(
  "Nuc ∩ eCLIP:",
  length(intersect(nuc_genes, common_eclip_genes)),
  "\n"
)

cat(
  "Total ∩ eCLIP:",
  length(intersect(total_genes, common_eclip_genes)),
  "\n"
)

cat(
  "Cyto ∩ Nuc ∩ Total:",
  length(
    Reduce(
      intersect,
      list(cyto_genes, nuc_genes, total_genes)
    )
  ),
  "\n"
)

cat(
  "Cyto ∩ Nuc ∩ Total ∩ eCLIP:",
  length(
    Reduce(
      intersect,
      list(
        cyto_genes,
        nuc_genes,
        total_genes,
        common_eclip_genes
      )
    )
  ),
  "\n"
)


# ============================================================
# Define Venn sets
# ============================================================

venn_sets <- list(
  Cyto = cyto_genes,
  Nuc = nuc_genes,
  Total = total_genes,
  eCLIP = common_eclip_genes
)


# ============================================================
# Calculate total number of unique genes in the union
# ============================================================

venn_union <- Reduce(
  union,
  venn_sets
)

n_union <- length(venn_union)

cat(
  "Total unique genes in union:",
  n_union,
  "\n"
)


# ============================================================
# Calculate percentages
#
# Percentage = number of genes in region / total genes
# in the union of all four sets × 100
# ============================================================

# Create labels showing:
#     gene count
#     percentage of total union
#
# VennDiagram will use these labels in the regions.

region_counts <- calculate.overlap(venn_sets)

region_labels <- sapply(
  region_counts,
  function(x) {
    
    n <- length(x)
    pct <- (n / n_union) * 100
    
    paste0(
      n,
      " (",
      sprintf("%.2f", pct),
      "%)"
    )
  }
)


# ============================================================
# Create Venn diagram
# ============================================================

venn_plot <- venn.diagram(
  
  x = venn_sets,
  filename = NULL,
  
  fill = c(
    "#08519C",
    "#800080",
    "#028A0F",
    "#BF6155"
  ),
  
  alpha = 0.5,
  
  col = "black",
  lwd = 1.2,
  
  # Display count + percentage
  print.mode = c("raw", "percent"),
  
  # Number formatting
  cex = 1.2,
  fontface = "bold",
  
  cat.cex = 1.2,
  cat.fontface = "bold",
  
  margin = 0.1,
  
  main = "Overlap of rMATS and Common eCLIP Genes",
  main.cex = 1.4,
  main.fontface = "bold"
)


# ============================================================
# Display Venn diagram
# ============================================================

grid.newpage()
grid.draw(venn_plot)


# ============================================================
# Save Venn diagram as PDF
# ============================================================

pdf(
  "Venn_rMATS_Cyto_Nuc_Total_common_eCLIP_genes.pdf",
  width = 8,
  height = 7
)

grid.newpage()
grid.draw(venn_plot)

dev.off()


# ============================================================
# Get genes common to all four sets
# ============================================================

all_four <- Reduce(
  intersect,
  list(
    cyto_genes,
    nuc_genes,
    total_genes,
    common_eclip_genes
  )
)


# ============================================================
# Print and save all-four overlap
# ============================================================

cat(
  "Genes common to Cyto + Nuc + Total + eCLIP:",
  length(all_four),
  "\n"
)

write.csv(
  data.frame(gene = all_four),
  "Cyto_Nuc_Total_eCLIP_common_genes.csv",
  row.names = FALSE
)


