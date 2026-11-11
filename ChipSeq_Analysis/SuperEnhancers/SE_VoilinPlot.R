#------------------------------------------
# SEs Violin Plot
#------------------------------------------
library(ggplot2)
se <- read.table("/home/umahesh284/Agastya_UM/M_ChipSeq/M_OP/batch2_t2p/DiffBind/summits_1000/K27ac/original_coordinates/oc/len/SE_MChip/diffbind_SE/Union/SE_union_named.bed", header=FALSE)
colnames(se) <- c("Chr","Start","End","SE_name")
nrow(se)
df <- read.csv("/home/umahesh284/Agastya_UM/M_ChipSeq/M_OP/batch2_t2p/DiffBind/summits_1000/K27ac/original_coordinates/oc/len/SE_MChip/diffbind_SE/Union/SE_significant.csv")
colnames(df) <- c("Chr","Start","End","width", "strand", "Conc", "Conc_HNRNPKsi", "Conc_Control", "Fold", "p.value", "FDR")

df <- merge(df, se,
            by=c("Chr","Start","End"),
            all.x=TRUE)
df

df$negLogFDR <- -log10(df$FDR)

df$significance <- "Not Significant"
df$significance[df$FDR < 0.05 & df$Fold >= 0.25]  <- "Enriched"
df$significance[df$FDR < 0.05 & df$Fold <= -0.25] <- "Depleted"

# -----------------------------
# Labels
# -----------------------------
df$label <- df$SE_name
df <- df[df$significance!="Not Significant",]
up <- df[df$significance=="Enriched",]
down <- df[df$significance=="Depleted",]
nrow(df)
# --------------------------
# Prepare data
# --------------------------

df <- data.frame(
  log2FC = df$Fold
)

df <- na.omit(df)
nrow(df)
# Sample size
n <- nrow(df)

label_text <- paste0(
  "n = ", format(n, big.mark = ",")
)

# Y-axis range
ymax <- max(df$log2FC, na.rm = TRUE)
ymin <- min(df$log2FC, na.rm = TRUE)
yrange <- ymax - ymin
# --------------------------
# Plot
# --------------------------

p <- ggplot(df, aes(x = "", y = log2FC)) +
  
  # Violin
  geom_violin(
    fill = "#D55E00",
    trim = FALSE,
    width = 0.8,
    alpha = 0.60,
    color = NA
  ) +
  
  # Boxplot
  geom_boxplot(
    width = 0.20,
    outlier.shape = NA,
    linewidth = 0.4,
    fill = "#D55E00",
    color = "black",
    alpha = 0.9
  ) +
  
  # Jitter points
  geom_jitter(
    width = 0.12,
    size = 0.35,
    alpha = 0.12,
    color = "black"
  ) +
  
  # Zero reference line
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.5
  ) +
  
  # Sample size
  annotate(
    "label",
    x = 0.65,
    y = ymax * 0.95,
    label = label_text,
    hjust = 1,
    vjust = 1,
    size = 4,
    fill = "white",
    label.size = 0.5
  ) +
  
  labs(
    x = NULL,
    y = expression(Log[2]~Fold~Change),
    title = "Distribution of Differential SEs Fold Changes"
  ) +
  
  theme_bw(base_size = 14) +
  
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    
    panel.border = element_rect(
      colour = "black",
      linewidth = 0.8
    ),
    
    axis.title.y = element_text(
      size = 16,
      face = "bold"
    ),
    
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    
    axis.text.y = element_text(
      size = 12,
      colour = "black"
    ),
    
    axis.ticks = element_line(
      linewidth = 0.8
    ),
    
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 16
    )
  ) +
  
  coord_cartesian(
    ylim = c(
      ymin - 0.5 * yrange,
      ymax + 0.5 * yrange
    ),
    clip = "on"
  )

p

# --------------------------
# Save
# --------------------------

ggsave(
  "/home/umahesh284/SEs_log2FC_publication.pdf",
  p,
  width = 4.5,
  height = 5,
  dpi = 600
)
