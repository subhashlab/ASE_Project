#-----------------------------------------------------------------------------------------
# *Create bigwig files with average signal across replicates for both H3K27ac and H3K27me3* 
#-----------------------------------------------------------------------------------------

#bigwigAverage -b /iitjfs/home/2025rbb1002/M_ChipSeq/batch2_raw/Peaks/H3K27ac/bedgraph/bw/Unknown_CR300-001R0001.cpm.bw /iitjfs/home/2025rbb1002/M_ChipSeq/batch2_raw/Peaks/H3K27ac/bedgraph/bw/Unknown_CR300-001R0002.cpm.bw -o ac_avg_cntrl.bw

#bigwigAverage -b /iitjfs/home/2025rbb1002/M_ChipSeq/batch2_raw/Peaks/H3K27ac/bedgraph/bw/Unknown_CR300-001R0003.cpm.bw /iitjfs/home/2025rbb1002/M_ChipSeq/batch2_raw/Peaks/H3K27ac/bedgraph/bw/Unknown_CR300-001R0004.cpm.bw /iitjfs/home/2025rbb1002/M_ChipSeq/batch2_raw/Peaks/H3K27ac/bedgraph/bw/Unknown_CR300-001R0005.cpm.bw /iitjfs/home/2025rbb1002/M_ChipSeq/batch2_raw/Peaks/H3K27ac/bedgraph/bw/Unknown_CR300-001R0006.cpm.bw -o ac_avg_hn.bw

bigwigAverage -b /iitjfs/home/2025rbb1002/M_ChipSeq/batch2_raw/Peaks/H3K27me3/bedgraph/bw/Unknown_CR300-001R0007.cpm.bw /iitjfs/home/2025rbb1002/M_ChipSeq/batch2_raw/Peaks/H3K27me3/bedgraph/bw/Unknown_CR300-001R0008.cpm.bw -o me3_avg_cntrl.bw

#bigwigAverage -b /iitjfs/home/2025rbb1002/M_ChipSeq/batch2_raw/Peaks/H3K27me3/bedgraph/bw/Unknown_CR300-001R0009.cpm.bw /iitjfs/home/2025rbb1002/M_ChipSeq/batch2_raw/Peaks/H3K27me3/bedgraph/bw/Unknown_CR300-001R0010.cpm.bw /iitjfs/home/2025rbb1002/M_ChipSeq/batch2_raw/Peaks/H3K27me3/bedgraph/bw/Unknown_CR300-001R0011.cpm.bw /iitjfs/home/2025rbb1002/M_ChipSeq/batch2_raw/Peaks/H3K27me3/bedgraph/bw
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------------------
# *Compute matrix for dsrip genes and TEs regions, cluster them into 2clusters AND plot heatmap for the cluster1 with the high magnitude of signal*
#------------------------------------------------------------------------------------------------------------------------------------------------
computeMatrix scale-regions   -S /iitjfs/home/2025rbb1002/M_ChipSeq/M_OP/batch2_t2p/DiffBind/summits_1000/integration/ASE/ASE_heatmaps/plots/ac_avg_cntrl.bw /iitjfs/home/2025rbb1002/M_ChipSeq/M_OP/batch2_t2p/DiffBind/summits_1000/integration/ASE/ASE_heatmaps/plots/ac_avg_hn.bw -b 2000 -a 2000 -R /iitjfs/home/2025rbb1002/M_ChipSeq/M_OP/batch2_t2p/DiffBind/summits_1000/integration/dsRIP_new/TEs/dsrip_TEs.bed --missingDataAsZero   -o /iitjfs/home/2025rbb1002/M_ChipSeq/M_OP/batch2_t2p/DiffBind/summits_1000/integration/dsRIP_new/TEs/dsrip_TEs_matrix.gz

plotHeatmap -m /iitjfs/home/2025rbb1002/M_ChipSeq/M_OP/batch2_t2p/DiffBind/summits_1000/integration/dsRIP_new/TEs/dsrip_TEs_matrix.gz --kmeans 2   --outFileSortedRegions /iitjfs/home/2025rbb1002/M_ChipSeq/M_OP/batch2_t2p/DiffBind/summits_1000/integration/dsRIP_new/TEs/clustered_regions_TEs_2.bed --colorMap RdBu_r   --samplesLabel   ac_Cntrl ac_hn -out /iitjfs/home/2025rbb1002/M_ChipSeq/M_OP/batch2_t2p/DiffBind/summits_1000/integration/dsRIP_new/TEs/heatmap_dsrip_TEs_2.pdf


computeMatrix scale-regions   -S /iitjfs/home/2025rbb1002/M_ChipSeq/M_OP/batch2_t2p/DiffBind/summits_1000/integration/ASE/ASE_heatmaps/plots/me3_avg_cntrl.bw /iitjfs/home/2025rbb1002/M_ChipSeq/M_OP/batch2_t2p/DiffBind/summits_1000/integration/ASE/ASE_heatmaps/plots/me3_avg_hn.bw -b 2000 -a 2000 -R /iitjfs/home/2025rbb1002/M_ChipSeq/M_OP/batch2_t2p/DiffBind/summits_1000/integration/dsRIP_new/genebody/withoutTT/withoutTT2/cluster1_withoutTT_new_me3.bed --missingDataAsZero -o /iitjfs/home/2025rbb1002/M_ChipSeq/M_OP/batch2_t2p/DiffBind/summits_1000/integration/dsRIP_new/genebody/withoutTT/withoutTT2/dsrip_genebody_me3_without_TT_cluster1_matrix.gz

plotHeatmap   -m /iitjfs/home/2025rbb1002/M_ChipSeq/M_OP/batch2_t2p/DiffBind/summits_1000/integration/dsRIP_new/genebody/withoutTT/withoutTT2/dsrip_genebody_me3_without_TT_cluster1_matrix.gz --kmeans 1   --outFileSortedRegions /iitjfs/home/2025rbb1002/M_ChipSeq/M_OP/batch2_t2p/DiffBind/summits_1000/integration/dsRIP_new/genebody/withoutTT/withoutTT2/clustered_regions_genebody_me3_withoutTT_cluster1.bed --colorMap RdBu_r   --samplesLabel   me3_Cntrl me3_hn -out /iitjfs/home/2025rbb1002/M_ChipSeq/M_OP/batch2_t2p/DiffBind/summits_1000/integration/dsRIP_new/genebody/withoutTT/withoutTT2/heatmap_genebody_me3_withoutTT_1.pdf
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#------------------------------------------------------------------
# *Selected genes with length above 5kb from dsrip enriched genes*
#------------------------------------------------------------------

dsrip <- read.csv("~/Agastya_UM/M_ChipSeq/M_OP/batch2_t2p/DiffBind/summits_1000/integration/dsRIP_new/dsrip_with_strand.bed",sep = "\t", header = F)
head(dsrip)
colnames(dsrip) <- c("chr","start","end","gene_id","pseduo_score","strand")
dsrip$length <- dsrip$end - dsrip$start
dsrip <- dsrip[dsrip$chr != "TT",]
dsrip
write.csv(dsrip,"~/Agastya_UM/M_ChipSeq/M_OP/batch2_t2p/DiffBind/summits_1000/integration/dsRIP_new/dsrip_without_TT_with_length.csv",row.names = F)
nrow(dsrip)

# Select only genes with length more than 5kb
dsrip_5kb <- dsrip[dsrip$length>5000,]
nrow(dsrip_5kb)
write.csv(dsrip_5kb,"~/Agastya_UM/M_ChipSeq/M_OP/batch2_t2p/DiffBind/summits_1000/integration/dsRIP_new/dsrip_genes_with_length_above_5kb.csv",row.names = F)
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#--------------------------------------------------------------------------------------------------------------------------
# *The >5kb length genes in cluster1 in both H3K27ac and H3K27me3 were used to get their 4 regions(Promoter,5'UTR,CDS,3'UTR)*
#--------------------------------------------------------------------------------------------------------------------------

library(GenomicFeatures)
library(GenomicRanges)
library(rtracklayer)
library(dplyr)

#--------------------------------------------------
# Read gene list
#--------------------------------------------------
genes <- read.csv("/iitjfs/home/2025rbb1002/M_ChipSeq/M_OP/batch2_t2p/DiffBind/summits_1000/integration/dsRIP_new/4regions/avg/dsrip_genes_with_length_above_5kb_cluster1_me3_2.csv")

gene_list <- unique(genes$gene_id)

#--------------------------------------------------
# Load GTF
#--------------------------------------------------
txdb <- makeTxDbFromGFF(
    "/iitjfs/home/2025rbb1002/hg38/Homo_sapiens.GRCh38.113.chr.gtf"
)

#--------------------------------------------------
# Longest transcript per gene
#--------------------------------------------------
tx <- transcripts(
    txdb,
    columns = c("tx_id","gene_id")
)

tx_df <- as.data.frame(tx)
tx_df$tx_width <- width(tx)

longest_tx <- tx_df %>%
    filter(gene_id %in% gene_list) %>%
    group_by(gene_id) %>%
    slice_max(tx_width, n = 1, with_ties = FALSE)

#--------------------------------------------------
# Transcript features
#--------------------------------------------------
utr5 <- fiveUTRsByTranscript(txdb)
utr3 <- threeUTRsByTranscript(txdb)
cds  <- cdsBy(txdb, by = "tx")

#--------------------------------------------------
# Keep transcripts having ALL regions
# (effectively protein-coding transcripts
# with annotated 5'UTR and 3'UTR)
#--------------------------------------------------
good_tx <- Reduce(
    intersect,
    list(
        as.character(longest_tx$tx_id),
        names(utr5),
        names(utr3),
        names(cds)
    )
)

# keep only those representative transcripts
longest_tx <- longest_tx[
    as.character(longest_tx$tx_id) %in% good_tx,
]

tx_long <- tx[
    tx$tx_id %in% longest_tx$tx_id
]

#--------------------------------------------------
# Report counts
#--------------------------------------------------
#--------------------------------------------------
# Report counts
#--------------------------------------------------
cat("Genes in input file:", length(gene_list), "\n")
cat("Protein-coding genes with 5'UTR, CDS and 3'UTR:",
    nrow(longest_tx), "\n")

# List final gene names
cat("\nFinal gene names:\n")
print(longest_tx$gene_id)
write.csv(longest_tx$gene_id,"final_pc_5k_genes_me3.csv")
#--------------------------------------------------
# Promoter (TSS-1kb -> TSS)
#--------------------------------------------------
promoter <- promoters(
    tx_long,
    upstream = 1000,
    downstream = 500
)

#--------------------------------------------------
# Collapse 5'UTR to one region per transcript
#--------------------------------------------------
utr5_region <- GRanges()

for(id in good_tx){

    x <- utr5[[id]]

    utr5_region <- c(
        utr5_region,
        GRanges(
            seqnames = seqnames(x)[1],
            ranges = IRanges(
                start = min(start(x)),
                end   = max(end(x))
            ),
            strand = strand(x)[1]
        )
    )
}

#--------------------------------------------------
# Collapse CDS to one span per transcript
# first CDS start -> last CDS end
#--------------------------------------------------
cds_region <- GRanges()

for(id in good_tx){

    x <- cds[[id]]

    cds_region <- c(
        cds_region,
        GRanges(
            seqnames = seqnames(x)[1],
            ranges = IRanges(
                start = min(start(x)),
                end   = max(end(x))
            ),
            strand = strand(x)[1]
        )
    )
}

#--------------------------------------------------
# Collapse 3'UTR to one region per transcript
#--------------------------------------------------
utr3_region <- GRanges()

for(id in good_tx){

    x <- utr3[[id]]

    utr3_region <- c(
        utr3_region,
        GRanges(
            seqnames = seqnames(x)[1],
            ranges = IRanges(
                start = min(start(x)),
                end   = max(end(x))
            ),
            strand = strand(x)[1]
        )
    )
}

#--------------------------------------------------
# Check counts
#--------------------------------------------------
cat("\nRegion counts:\n")
cat("Promoter :", length(promoter), "\n")
cat("5UTR     :", length(utr5_region), "\n")
cat("CDS span :", length(cds_region), "\n")
cat("3UTR     :", length(utr3_region), "\n")

#--------------------------------------------------
# Export BED files
#--------------------------------------------------
export(promoter,    "Promoter_TSSminus1kb_new_me3.bed")
export(utr5_region, "UTR5_region_new_me3.bed")
export(cds_region,  "CDS_span_new_me3.bed")
export(utr3_region, "UTR3_region_new_me3.bed")
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# *The distribution of genes and TEs in cluster1 were ploted in proportion plots*
#------------------------------------------------------------------------------

library(rtracklayer)
library(dplyr)
library(ggplot2)

# =====================================================
# LOAD GTF
# =====================================================

gtf <- import("/iitjfs/home/2025rbb1002/hg38/Homo_sapiens.GRCh38.113.chr.gtf")

genes <- gtf[gtf$type == "gene"]

# detect correct biotype column
biotype_col <- if ("gene_biotype" %in% colnames(mcols(genes))) {
  "gene_biotype"
} else {
  "gene_type"
}

# annotation table
gene_annot <- data.frame(
  gene_id = sub("\\..*", "", genes$gene_id),
  biotype = mcols(genes)[[biotype_col]]
)

gene_annot <- distinct(gene_annot, gene_id, .keep_all = TRUE)

# =====================================================
# FILES + COLUMN NAMES
# =====================================================

files <- c(
  "/iitjfs/home/2025rbb1002/M_ChipSeq/M_OP/batch2_t2p/DiffBind/summits_1000/integration/dsRIP_new/4regions/avg/dsrip_genes_with_length_above_5kb_cluster1_me3_2.csv"
  #"/iitjfs/home/2025rbb1002/M_ChipSeq/M_OP/batch2_t2p/DiffBind/summits_1000/K27ac/original_coordinates/oc/5k/Annotated_diff_peaks_K27ac_sig_gene_up_5kb_allgenes_detailed.csv"
)

# names for plot
# names for plot (based on unique genes)

sample_names <- c(
  paste0(
    "Depleted (",
    length(unique(sub("\\..*", "", read.csv(files[1], stringsAsFactors = FALSE)[["gene_id"]]))),
    " genes)"
  ),
  paste0(
    "Enriched (",
    length(unique(sub("\\..*", "", read.csv(files[2], stringsAsFactors = FALSE)[["gene_id"]]))),
    " genes)"
  )
)

sample_names <- c(
  paste0(
    "H3K27ac (",
    length(unique(sub("\\..*", "", read.csv(files[1], stringsAsFactors = FALSE)[["gene_id"]]))),
    " genes)"
  )
)

# column containing gene IDs
gene_column <- "gene_id"

# =====================================================
# PROCESS ALL FILES
# =====================================================

all_plot_df <- list()

for (i in seq_along(files)) {
  
  f <- files[i]
  
  # -----------------------------
  # READ FILE
  # -----------------------------
  
  df <- read.csv(f, stringsAsFactors = FALSE)
  
  # -----------------------------
  # EXTRACT GENE IDS
  # -----------------------------
  
  gene_list <- unique(sub("\\..*", "", df[[gene_column]]))
  
  # -----------------------------
  # MAP BIOTYPE
  # -----------------------------
  
  mapped <- data.frame(gene_id = gene_list) %>%
    left_join(gene_annot, by = "gene_id")
  
  mapped$biotype[is.na(mapped$biotype)] <- "Unknown"
  
  # -----------------------------
  # COUNT + PROPORTION
  # -----------------------------
  
  plot_df <- mapped %>%
    group_by(biotype) %>%
    summarise(count = n(), .groups = "drop") %>%
    mutate(freq = count / sum(count))
  
  # custom sample name
  plot_df$sample <- sample_names[i]
  
  all_plot_df[[i]] <- plot_df
}

# =====================================================
# COMBINE
# =====================================================

final_plot_df <- bind_rows(all_plot_df)

# keep order
final_plot_df$sample <- factor(
  final_plot_df$sample,
  levels = sample_names
)

# =====================================================
# PLOT
# =====================================================

p <- ggplot(final_plot_df,
            aes(x = sample,
                y = freq,
                fill = biotype)) +
  
  geom_bar(stat = "identity") +
  
  scale_fill_manual(
    name = "Biotype",
    values = c(
      lncRNA                           = "#4E79A7",
      miRNA                            = "#E15759",
      misc_RNA                         = "#76B7B2",
      processed_pseudogene             = "#9C755F",
      protein_coding                   = "#F28E2B",
      rRNA_pseudogene                  = "#EDC948",
      snRNA                            = "#B07AA1",
      TEC                              = "#FF9DA7",
      transcribed_processed_pseudogene = "#59A14F",
      transcribed_unitary_pseudogene   = "#8CD17D",
      unitary_pseudogene               = "#BAB0AC",
      unprocessed_pseudogene           = "#D37295",
      Unknown                          = "#999999"
    )
  ) +
  
  ylab("Proportion") +
  xlab("") +
  
  theme_minimal() +
  
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5)
  )

print(p)

# =====================================================
# SAVE
# =====================================================

dir.create("plots", showWarnings = FALSE)

ggsave("/iitjfs/home/2025rbb1002/M_ChipSeq/M_OP/batch2_t2p/DiffBind/summits_1000/integration/dsRIP_new/4regions/avg/biotype_stacked_K27me3_rip_chip_110.pdf",
       plot = p,
       width = 3,
       height = 5)

write.csv(mapped,"/iitjfs/home/2025rbb1002/M_ChipSeq/M_OP/batch2_t2p/DiffBind/summits_1000/integration/dsRIP_new/4regions/avg/H3K27me3_gene_biotypes_110.csv")   
#####################################
# TEs
####################################
library(rtracklayer)
library(dplyr)
library(ggplot2)

gtf <- import("/iitjfs/home/2025rbb1002/hg38/GRCh38_Ensembl_rmsk_TE.gtf")

genes <- gtf[gtf$type == "exon"]

# detect correct biotype column
biotype_col <- if ("class_id" %in% colnames(mcols(genes))) {
  "class_id"
} else {
  "class_type"
}

# annotation table
gene_annot <- data.frame(
  transcript_id = sub("\\..*", "", genes$transcript_id),
  class_id = mcols(genes)[[biotype_col]]
)

gene_annot <- distinct(gene_annot, transcript_id, .keep_all = TRUE)
file <- read.csv("/iitjfs/home/2025rbb1002/M_ChipSeq/M_OP/batch2_t2p/DiffBind/summits_1000/integration/dsRIP_new/TEs/clustered_regions_TEs_me3_1of2.bed",sep="\t")
colnames(file) <- c("chrom","start","end","transcript_id","score","strand","thickStart","thickEnd","itemRGB", "blockCount", "blockSizes", "blockStart", "deepTools_group")
file_annotated <- left_join(file,gene_annot, by = "transcript_id")
write.csv(file_annotated,"/iitjfs/home/2025rbb1002/M_ChipSeq/M_OP/batch2_t2p/DiffBind/summits_1000/integration/dsRIP_new/TEs/clustered_regions_TEs_me3_1of2_TE_class_annotated.csv")
##########################################
# Plots
##########################################
library(dplyr)
library(ggplot2)

# =========================================
# LOAD FILES
# =========================================
ripNac_down <- read.csv(
  "/iitjfs/home/2025rbb1002/M_ChipSeq/M_OP/batch2_t2p/DiffBind/summits_1000/integration/dsRIP_new/TEs/clustered_regions_TEs_1of2_TE_class_annotated.csv"
)

# =========================================
# SIMPLIFY CLASS IDS
# =========================================

simplify_class <- function(x) {
  
  case_when(
    grepl("^DNA", x)        ~ "DNA",
    grepl("^LTR", x)        ~ "LTR",
    grepl("^LINE", x)       ~ "LINE",
    grepl("^SINE", x)       ~ "SINE",
    grepl("^RC\\?", x)      ~ "RC?",
    grepl("^Retroposon", x) ~ "Retroposon",
    grepl("^Satellite", x)  ~ "Satellite",
    grepl("^Unknown", x)    ~ "Unknown",
    TRUE ~ "Other"
  )
}
ripNac_down <- ripNac_down %>%
  distinct(transcript_id, .keep_all = TRUE)

ripNac_down$class_group <- simplify_class(ripNac_down$class_id)

# =========================================
# REMOVE "Other"
# =========================================

ripNac_down <- ripNac_down %>%
  filter(class_group != "Other")

# =========================================
# COUNT PROPORTIONS
# =========================================

p1 <- ripNac_down %>%
  group_by(class_group) %>%
  summarise(count = n(), .groups = "drop") %>%
  mutate(
    freq = count / sum(count),
    sample = paste0("Depleted(", nrow(ripNac_down), " TEs)")
  )

# =========================================
# COMBINE
# =========================================

plot_df <- bind_rows(p1)

# =========================================
# ORDER GROUPS
# =========================================

plot_df$class_group <- factor(
  plot_df$class_group,
  levels = c(
    "DNA",
    "LINE",
    "LTR",
    "RC?",
    "Retroposon",
    "Satellite",
    "SINE",
    "Unknown"
  )
)

# =========================================
# PLOT
# =========================================

p <- ggplot(
  plot_df,
  aes(
    x = sample,
    y = freq,
    fill = class_group
  )
) +
  
  geom_bar(stat = "identity") +
  
  scale_fill_manual(
    name = "class",
    values = c(
      DNA        = "#3893C4",  # deep raspberry
      LINE       = "#A32C5A",  # dusty pink
      LTR        = "#0A8068",  # deep teal-green
      `RC?`      = "#5D79BC",  # muted blue
      Retroposon = "#BAD988",  # light green
      Satellite  = "#8E54A1",  # purple
      SINE       = "#E8A4B2",  # blue-cyan
      Unknown    = "#B0B0B0"   # gray
    )
  ) +
  
  ylab("Proportion") +
  xlab("") +
  
  theme_minimal() +
  
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5)
  )

print(p)

# =========================================
# SAVE
# =========================================

ggsave(
  "/iitjfs/home/2025rbb1002/M_ChipSeq/M_OP/batch2_t2p/DiffBind/summits_1000/integration/dsRIP_new/TEs/clustered_regions_TEs_ac_1of2_TE_class.pdf",
  plot = p,
  width = 6,
  height = 5,
  units = "in"
)

#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#-----------------------------------
# *Mean signal over all four regions*
#-----------------------------------

multiBigwigSummary BED-file \
  --bwfiles /iitjfs/home/2025rbb1002/M_ChipSeq/M_OP/batch2_t2p/DiffBind/summits_1000/integration/ASE/ASE_heatmaps/plots/me3_avg_cntrl.bw /iitjfs/home/2025rbb1002/M_ChipSeq/M_OP/batch2_t2p/DiffBind/summits_1000/integration/ASE/ASE_heatmaps/plots/me3_avg_hn.bw \
  --BED /iitjfs/home/2025rbb1002/M_ChipSeq/M_OP/batch2_t2p/DiffBind/summits_1000/integration/dsRIP_new/genebody/withoutTT/withoutTT2/Promoter_TSSminus1kb_new_me3.bed \
  --outRawCounts /iitjfs/home/2025rbb1002/M_ChipSeq/M_OP/batch2_t2p/DiffBind/summits_1000/integration/dsRIP_new/genebody/withoutTT/withoutTT2/Promoter_signal_me3_new2.tab \
  -out /iitjfs/home/2025rbb1002/M_ChipSeq/M_OP/batch2_t2p/DiffBind/summits_1000/integration/dsRIP_new/genebody/withoutTT/withoutTT2/Promoter_me3_new2.npz
  
#UTR5 

multiBigwigSummary BED-file \
  --bwfiles /iitjfs/home/2025rbb1002/M_ChipSeq/M_OP/batch2_t2p/DiffBind/summits_1000/integration/ASE/ASE_heatmaps/plots/me3_avg_cntrl.bw /iitjfs/home/2025rbb1002/M_ChipSeq/M_OP/batch2_t2p/DiffBind/summits_1000/integration/ASE/ASE_heatmaps/plots/me3_avg_hn.bw \
  --BED /iitjfs/home/2025rbb1002/M_ChipSeq/M_OP/batch2_t2p/DiffBind/summits_1000/integration/dsRIP_new/genebody/withoutTT/withoutTT2/UTR5_region_new_me3.bed \
  --outRawCounts /iitjfs/home/2025rbb1002/M_ChipSeq/M_OP/batch2_t2p/DiffBind/summits_1000/integration/dsRIP_new/genebody/withoutTT/withoutTT2/UTR5_region_new_me3.tab \
  -out /iitjfs/home/2025rbb1002/M_ChipSeq/M_OP/batch2_t2p/DiffBind/summits_1000/integration/dsRIP_new/genebody/withoutTT/withoutTT2/UTR5_region_new_me3.npz
  
#CDS
multiBigwigSummary BED-file \
  --bwfiles /iitjfs/home/2025rbb1002/M_ChipSeq/M_OP/batch2_t2p/DiffBind/summits_1000/integration/ASE/ASE_heatmaps/plots/me3_avg_cntrl.bw /iitjfs/home/2025rbb1002/M_ChipSeq/M_OP/batch2_t2p/DiffBind/summits_1000/integration/ASE/ASE_heatmaps/plots/me3_avg_hn.bw \
  --BED /iitjfs/home/2025rbb1002/M_ChipSeq/M_OP/batch2_t2p/DiffBind/summits_1000/integration/dsRIP_new/genebody/withoutTT/withoutTT2/CDS_span_new_me3.bed \
  --outRawCounts /iitjfs/home/2025rbb1002/M_ChipSeq/M_OP/batch2_t2p/DiffBind/summits_1000/integration/dsRIP_new/genebody/withoutTT/withoutTT2/CDS_span_new_me3.tab \
  -out /iitjfs/home/2025rbb1002/M_ChipSeq/M_OP/batch2_t2p/DiffBind/summits_1000/integration/dsRIP_new/genebody/withoutTT/withoutTT2/CDS_span_new_me3.npz
  
  
#UTR3
multiBigwigSummary BED-file \
  --bwfiles /iitjfs/home/2025rbb1002/M_ChipSeq/M_OP/batch2_t2p/DiffBind/summits_1000/integration/ASE/ASE_heatmaps/plots/me3_avg_cntrl.bw /iitjfs/home/2025rbb1002/M_ChipSeq/M_OP/batch2_t2p/DiffBind/summits_1000/integration/ASE/ASE_heatmaps/plots/me3_avg_hn.bw \
  --BED /iitjfs/home/2025rbb1002/M_ChipSeq/M_OP/batch2_t2p/DiffBind/summits_1000/integration/dsRIP_new/genebody/withoutTT/withoutTT2/UTR3_region_new_me3.bed \
  --outRawCounts /iitjfs/home/2025rbb1002/M_ChipSeq/M_OP/batch2_t2p/DiffBind/summits_1000/integration/dsRIP_new/genebody/withoutTT/withoutTT2/UTR3_region_new_me3.tab \
  -out /iitjfs/home/2025rbb1002/M_ChipSeq/M_OP/batch2_t2p/DiffBind/summits_1000/integration/dsRIP_new/genebody/withoutTT/withoutTT2/UTR3_region_new_me3.npz
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#--------------------------------------------------------------
# *Statistical analysis of mean signal difference in each region*
#--------------------------------------------------------------

library(tidyr)
library(dplyr)
library(ggplot2)

#---------------------------------------------------------
# Function
#---------------------------------------------------------

read_region <- function(file, region){

  x <- read.table(
    file,
    header = FALSE,
    colClasses = c(
      "character",
      "character",
      "character",
      "numeric",
      "numeric"
    )
  )

  colnames(x)[4:5] <- c("Control", "Treatment")
#0.3514772
#0.05513738

#0.04782976_ac_up
#0.1555187_ac_down
  x$Control   <- log2(x$Control + 0.1555187)
  x$Treatment <- log2(x$Treatment + 0.1555187)

  # Mean difference
  mean.diff <- mean(x$Treatment - x$Control)

  # One-sided paired Wilcoxon
  wt <- wilcox.test(
    x$Treatment,
    x$Control,
    paired = TRUE,
    alternative = "greater"
  )

  ymax <- max(c(x$Control, x$Treatment))

  # Long format
  long <- pivot_longer(
    x,
    cols = c(Control, Treatment),
    names_to = "Group",
    values_to = "Signal"
  )

  long$Region <- region

  cat("\n==============================\n")
  cat(region, "\n")
  cat("==============================\n")
  cat("Mean Control   :", mean(x$Control), "\n")
  cat("Mean Treatment :", mean(x$Treatment), "\n")
  cat("Mean Difference:", mean.diff, "\n")
  print(wt)

  list(
    data = long,
    stats = data.frame(
      Region = region,
      Mean_Control = mean(x$Control),
      Mean_Treatment = mean(x$Treatment),
      Mean_Difference = mean.diff,
      pvalue = wt$p.value,
      ymax = ymax
    )
  )
}

#---------------------------------------------------------
# Read all regions
#---------------------------------------------------------

prom <- read_region(
  "Promoter_signal_me3_new2.tab",
  "PROMOTER"
)

utr5 <- read_region(
  "UTR5_region_new_me3.tab",
  "UTR5"
)

cds <- read_region(
  "CDS_region_new_me3.tab",
  "CDS"
)

utr3 <- read_region(
  "UTR3_region_new_me3.tab",
  "UTR3"
)

#---------------------------------------------------------
# Combine
#---------------------------------------------------------

df <- bind_rows(
  prom$data,
  utr5$data,
  cds$data,
  utr3$data
)

stats <- bind_rows(
  prom$stats,
  utr5$stats,
  cds$stats,
  utr3$stats
)

df$Region <- factor(
  df$Region,
  levels = c("PROMOTER","UTR5","CDS","UTR3")
)

stats$Region <- factor(
  stats$Region,
  levels = c("PROMOTER","UTR5","CDS","UTR3")
)

stats$label <- paste0(
  "mean difference = ",
  round(stats$Mean_Difference,3),
  "\n",
  "p = ",
  signif(stats$pvalue,3)
)

#---------------------------------------------------------
# Plot
#---------------------------------------------------------

p <- ggplot(
  df,
  aes(Group, Signal, fill = Group)
) +

  geom_violin(
    trim = FALSE,
    alpha = 0.6
  ) +

  geom_boxplot(
    width = 0.12,
    fill = "white",
    outlier.shape = NA
  ) +

  facet_wrap(
    ~Region,
    nrow = 1
  ) +

  geom_text(
    data = stats,
    aes(
      x = 1.5,
      y = ymax * 1.05,
      label = label
    ),
    inherit.aes = FALSE,
    size = 3
  ) +

  theme_bw() +

  theme(
    legend.position = "none",
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5)
  ) +

  labs(
    title = "Average H3K27me3 signal across dsRIP-enriched genes",
    x = "",
    y = "log2(CPM)"
  )

print(p)

#---------------------------------------------------------
# Save
#---------------------------------------------------------

ggsave(
  "Average_H3K27me3_signal_regions_mean_genes.pdf",
  plot = p,
  width = 7,
  height = 4
)

#---------------------------------------------------------
# Summary table
#---------------------------------------------------------

print(stats)

write.csv(
  stats,
  "Average_H3K27me3_Wilcoxon_Results_mean_genes.csv",
  row.names = FALSE
)

#---------------------------------------
# minimum value to add to avoid infinity
#---------------------------------------
promoter <- read.csv("~/Agastya_UM/M_ChipSeq/M_OP/batch2_t2p/DiffBind/summits_1000/integration/dsRIP_new/genebody/withoutMT/withoutMT2/diffbind/24_4regions/Promoter_signal_me3_new2",sep="\t",header = T)
utr3 <- read.csv("~/Agastya_UM/M_ChipSeq/M_OP/batch2_t2p/DiffBind/summits_1000/integration/dsRIP_new/genebody/withoutMT/withoutMT2/diffbind/24_4regions/UTR3_region_new_me3.tab",sep="\t",header = T)
utr5 <- read.csv("~/Agastya_UM/M_ChipSeq/M_OP/batch2_t2p/DiffBind/summits_1000/integration/dsRIP_new/genebody/withoutMT/withoutMT2/diffbind/24_4regions/UTR5_region_new_me3.tab",sep = "\t",header = T)
cds <- read.csv("~/Agastya_UM/M_ChipSeq/M_OP/batch2_t2p/DiffBind/summits_1000/integration/dsRIP_new/genebody/withoutMT/withoutMT2/diffbind/24_4regions/CDS_region_new_me3.tab",sep = "\t",header = T)

promoter
colnames(promoter) <- c("chr","start","end","cntrl","hn")
colnames(utr3) <- c("chr","start","end","cntrl","hn")
colnames(utr5) <- c("chr","start","end","cntrl","hn")
colnames(cds) <- c("chr","start","end","cntrl","hn")
min(promoter$cntrl)
min(promoter$hn)
mins <- c(min(promoter$cntrl),min(promoter$hn),min(utr3$cntrl),min(utr3$hn),min(utr5$cntrl),min(utr5$hn),min(cds$cntrl),min(cds$hn))
mins
min(mins)
sort(mins)

DNA <- read.csv("~/Agastya_UM/M_ChipSeq/M_OP/batch2_t2p/DiffBind/summits_1000/integration/dsRIP_new/TEs/DNA_me3.tab",sep="\t",header = T)
LTR <- read.csv("~/Agastya_UM/M_ChipSeq/M_OP/batch2_t2p/DiffBind/summits_1000/integration/dsRIP_new/TEs/LTR_me3.tab",sep="\t",header = T)
LINE <- read.csv("~/Agastya_UM/M_ChipSeq/M_OP/batch2_t2p/DiffBind/summits_1000/integration/dsRIP_new/TEs/LINE_me3.tab",sep = "\t",header = T)
SINE <- read.csv("~/Agastya_UM/M_ChipSeq/M_OP/batch2_t2p/DiffBind/summits_1000/integration/dsRIP_new/TEs/SINE_me3.tab",sep = "\t",header = T)

DNA
colnames(DNA) <- c("chr","start","end","cntrl","hn")
colnames(LTR) <- c("chr","start","end","cntrl","hn")
colnames(LINE) <- c("chr","start","end","cntrl","hn")
colnames(SINE) <- c("chr","start","end","cntrl","hn")
min(DNA$cntrl)
min(DNA$hn)
mins <- c(min(DNA$cntrl),min(DNA$hn),min(LTR$cntrl),min(LTR$hn),min(LINE$cntrl),min(LINE$hn),min(SINE$cntrl),min(SINE$hn))
mins
min(mins)
sort(mins)

#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#----------------------------------------
# *Separating the TEs based on their class*
#----------------------------------------

dsrip <- read.csv("~/Agastya_UM/M_ChipSeq/M_OP/batch2_t2p/DiffBind/summits_1000/integration/dsRIP_new/TEs/dsrip_TEs.bed", sep = "\t", header = F)
cluster1 <- read.csv("~/Agastya_UM/M_ChipSeq/M_OP/batch2_t2p/DiffBind/summits_1000/integration/dsRIP_new/TEs/dsrip_TEs_cluster1_me3.bed", header = F,sep = "\t")
cluster1
dsrip

library(dplyr)
cluster1_dsrip <- inner_join(dsrip,cluster1, by = "V4")
factor(cluster1_dsrip$V9.x)

#DNA LINE LTR Satellite SINE
DNA <- cluster1_dsrip[cluster1_dsrip$V9.x=="DNA",]
LINE <- cluster1_dsrip[cluster1_dsrip$V9.x=="LINE",]
LTR <- cluster1_dsrip[cluster1_dsrip$V9.x=="LTR",]
Satellite <- cluster1_dsrip[cluster1_dsrip$V9.x=="Satellite",]
SINE <- cluster1_dsrip[cluster1_dsrip$V9.x=="SINE",]
Retroposon <- cluster1_dsrip[cluster1_dsrip$V9.x=="Retroposon",]

write.table(DNA,"~/Agastya_UM/M_ChipSeq/M_OP/batch2_t2p/DiffBind/summits_1000/integration/dsRIP_new/TEs/DNA_me3.bed")
write.table(LINE,"~/Agastya_UM/M_ChipSeq/M_OP/batch2_t2p/DiffBind/summits_1000/integration/dsRIP_new/TEs/LINE_me3.bed")
write.table(LTR,"~/Agastya_UM/M_ChipSeq/M_OP/batch2_t2p/DiffBind/summits_1000/integration/dsRIP_new/TEs/LTR_me3.bed")
write.table(Satellite,"~/Agastya_UM/M_ChipSeq/M_OP/batch2_t2p/DiffBind/summits_1000/integration/dsRIP_new/TEs/Satellite_me3.bed")
write.table(SINE,"~/Agastya_UM/M_ChipSeq/M_OP/batch2_t2p/DiffBind/summits_1000/integration/dsRIP_new/TEs/SINE_me3.bed")
write.table(Retroposon,"~/Agastya_UM/M_ChipSeq/M_OP/batch2_t2p/DiffBind/summits_1000/integration/dsRIP_new/TEs/Retroposon_me3.bed")

#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------------
# *Checking the genes that have the enriched H3K27ac mean signal and depleted H3K27me3 signal*
#-------------------------------------------------------------------------------------------
################################################################################
# all dsrip genes
################################################################################
ac <- read.csv("~/Agastya_UM/M_ChipSeq/M_OP/batch2_t2p/DiffBind/summits_1000/integration/dsRIP_new/dsrip_withoutTT_ac.tab",sep="\t",header = T)
me3 <- read.csv("~/Agastya_UM/M_ChipSeq/M_OP/batch2_t2p/DiffBind/summits_1000/integration/dsRIP_new/dsrip_withoutTT_me3.tab",sep="\t",header = T)

ac
colnames(ac) <- c("chr","start","end","cntrl","hn")
colnames(me3) <- c("chr","start","end","cntrl","hn")
ac$diff <- ac$hn-ac$cntrl
me3$diff <- me3$hn-me3$cntrl
me3
ac
nrow(ac[ac$diff>0,])
nrow(ac[ac$diff<0,])
nrow(me3[me3$diff>0,])
nrow(me3[me3$diff<0,])

dsrip <- read.csv("~/Agastya_UM/M_ChipSeq/M_OP/batch2_t2p/DiffBind/summits_1000/integration/dsRIP_new/dsrip_withoutTT.bed",sep="\t",header = F)
dsrip
colnames(dsrip) <- c("chr","start","end","name","score","strand")
dsrip

ac_gene <- merge(ac,dsrip,by=c("chr","start","end"))
me3_gene <- merge(me3,dsrip,by=c("chr","start","end"))

ac_up <- ac_gene[ac_gene$diff>0,]
ac_down <- ac_gene[ac_gene$diff<0,]

me3_up <- me3_gene[me3_gene$diff>0,]
me3_down <- me3_gene[me3_gene$diff<0,]

up_up <- merge(ac_up,me3_up,by="name")
up_down <- merge(ac_up,me3_down,by="name")
down_up <- merge(ac_down,me3_up,by="name")
down_down <- merge(ac_down,me3_down,by="name")

nrow(up_up)
nrow(up_down)
nrow(down_up)
nrow(down_down)
write.csv(ac_up,"~/Agastya_UM/M_ChipSeq/M_OP/batch2_t2p/DiffBind/summits_1000/integration/dsRIP_new/genebody/withoutTT/withoutTT2/ac_up_dsrip.csv")
write.csv(ac_down,"~/Agastya_UM/M_ChipSeq/M_OP/batch2_t2p/DiffBind/summits_1000/integration/dsRIP_new/genebody/withoutTT/withoutTT2/ac_down_dsrip.csv")
write.csv(me3_up,"~/Agastya_UM/M_ChipSeq/M_OP/batch2_t2p/DiffBind/summits_1000/integration/dsRIP_new/genebody/withoutTT/withoutTT2/me3_up_dsrip.csv")
write.csv(me3_down,"~/Agastya_UM/M_ChipSeq/M_OP/batch2_t2p/DiffBind/summits_1000/integration/dsRIP_new/genebody/withoutTT/withoutTT2/me3_down_dsrip.csv")

write.csv(up_up,"~/Agastya_UM/M_ChipSeq/M_OP/batch2_t2p/DiffBind/summits_1000/integration/dsRIP_new/genebody/withoutTT/withoutTT2/ac_me3_up_up_dsrip.csv")
write.csv(up_down,"~/Agastya_UM/M_ChipSeq/M_OP/batch2_t2p/DiffBind/summits_1000/integration/dsRIP_new/genebody/withoutTT/withoutTT2/ac_me3_up_down_dsrip.csv")
write.csv(down_up,"~/Agastya_UM/M_ChipSeq/M_OP/batch2_t2p/DiffBind/summits_1000/integration/dsRIP_new/genebody/withoutTT/withoutTT2/ac_me3_down_up_dsrip.csv")
write.csv(down_down,"~/Agastya_UM/M_ChipSeq/M_OP/batch2_t2p/DiffBind/summits_1000/integration/dsRIP_new/genebody/withoutTT/withoutTT2/ac_me3_down_down_dsrip.csv")

################################
# TEs all
###############################
ac <- read.csv("~/Agastya_UM/M_ChipSeq/M_OP/batch2_t2p/DiffBind/summits_1000/integration/dsRIP_new/dsrip_TEs_ac.tab",sep="\t",header = T)
me3 <- read.csv("~/Agastya_UM/M_ChipSeq/M_OP/batch2_t2p/DiffBind/summits_1000/integration/dsRIP_new/dsrip_TEs_me3.tab",sep="\t",header = T)

nrow(me3)
colnames(ac) <- c("chr","start","end","cntrl","hn")
colnames(me3) <- c("chr","start","end","cntrl","hn")
ac$diff <- ac$hn-ac$cntrl
me3$diff <- me3$hn-me3$cntrl
me3
ac
nrow(ac[ac$diff==0,])
nrow(ac[ac$diff<0,])
nrow(me3[me3$diff==0,])
nrow(me3[me3$diff<0,])

dsrip <- read.csv("~/Agastya_UM/M_ChipSeq/M_OP/batch2_t2p/DiffBind/summits_1000/integration/dsRIP_new/TEs/dsrip_TEs.bed",sep="\t",header = F)
dsrip
colnames(dsrip) <- c("chr","start","end","name","score","strand","family","subfamily","class")
dsrip

ac_gene <- merge(ac,dsrip,by=c("chr","start","end"))
me3_gene <- merge(me3,dsrip,by=c("chr","start","end"))

ac_up <- ac_gene[ac_gene$diff>0,]
max(ac_up$diff)
ac_down <- ac_gene[ac_gene$diff<0,]
min(ac_down$diff)
me3_up <- me3_gene[me3_gene$diff>0,]
me3_up <- na.omit(me3_up)
max(me3_up$diff)
me3_down <- me3_gene[me3_gene$diff<0,]
me3_down <- na.omit(me3_down)
min(me3_down$diff)

sum(ac_up$diff)
sum(ac_down$diff)

sum(me3_up$diff)
sum(me3_down$diff)

up_up <- merge(ac_up,me3_up,by="name")
up_down <- merge(ac_up,me3_down,by="name")
down_up <- merge(ac_down,me3_up,by="name")
down_down <- merge(ac_down,me3_down,by="name")

nrow(up_up)
nrow(up_down)
nrow(down_up)
nrow(down_down)

write.csv(ac_up,"~/Agastya_UM/M_ChipSeq/M_OP/batch2_t2p/DiffBind/summits_1000/integration/dsRIP_new/genebody/withoutTT/withoutTT2/ac_up_TEs.csv")
write.csv(ac_down,"~/Agastya_UM/M_ChipSeq/M_OP/batch2_t2p/DiffBind/summits_1000/integration/dsRIP_new/genebody/withoutTT/withoutTT2/ac_down_TEs.csv")
write.csv(me3_up,"~/Agastya_UM/M_ChipSeq/M_OP/batch2_t2p/DiffBind/summits_1000/integration/dsRIP_new/genebody/withoutTT/withoutTT2/me3_up_TEs.csv")
write.csv(me3_down,"~/Agastya_UM/M_ChipSeq/M_OP/batch2_t2p/DiffBind/summits_1000/integration/dsRIP_new/genebody/withoutTT/withoutTT2/me3_down_TEs.csv")

write.csv(up_up,"~/Agastya_UM/M_ChipSeq/M_OP/batch2_t2p/DiffBind/summits_1000/integration/dsRIP_new/genebody/withoutTT/withoutTT2/ac_me3_up_up_TEs.csv")
write.csv(up_down,"~/Agastya_UM/M_ChipSeq/M_OP/batch2_t2p/DiffBind/summits_1000/integration/dsRIP_new/genebody/withoutTT/withoutTT2/ac_me3_up_down_TEs.csv")
write.csv(down_up,"~/Agastya_UM/M_ChipSeq/M_OP/batch2_t2p/DiffBind/summits_1000/integration/dsRIP_new/genebody/withoutTT/withoutTT2/ac_me3_down_up_TEs.csv")
write.csv(down_down,"~/Agastya_UM/M_ChipSeq/M_OP/batch2_t2p/DiffBind/summits_1000/integration/dsRIP_new/genebody/withoutTT/withoutTT2/ac_me3_down_down_TEs.csv")
################################################################################
# Proportion plots of genes and TEs with and without pattern
################################################################################

ac <- read.csv("~/Agastya_UM/M_ChipSeq/M_OP/batch2_t2p/DiffBind/summits_1000/integration/dsRIP_new/genebody/withoutTT/withoutTT2/dsrip_genes_propotions.tsv",sep="\t",row.names = 1)
prop_ac <- prop.table(as.matrix(ac), margin = 2)
ac

library(ggplot2)
ggplot(as.data.frame(as.table(prop_ac)),
       aes(x = Var2, y = Freq, fill = Var1)) +
  geom_col(width = 0.7) +
  scale_fill_manual(values = c(
    "High-Signal" = "#BA0900",
    "Low-Signal" = "#BA8E23"
  )) +
  theme_classic()

install.packages("VennDiagram")
library(VennDiagram)
library(grid)

# Enter numbers
A <- 1560
B <- 1610
AB <- 701

# Create 2-set Venn diagram
venn.plot <- draw.pairwise.venn(
  area1 = A,
  area2 = B,
  cross.area = AB,
  category = c("A", "B"),
  fill = c("#BA8E23", "#BA8E23"),
  alpha = 0.4,
  lwd = 2,
  cex = 1.5,
  cat.cex = 1.5,
  scaled = FALSE
)

grid.newpage()
grid.draw(venn.plot)
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
