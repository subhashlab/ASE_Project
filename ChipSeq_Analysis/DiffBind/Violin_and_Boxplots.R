#----------------------------------
# Voilin plot of significant peaks
#-----------------------------------
library(ggplot2)
library(dplyr)

# --------------------------
# Read data
# --------------------------

k27ac_up   <- read.csv("~/Agastya_UM/M_ChipSeq/M_OP/batch2_t2p/DiffBind/summits_1000/K27ac/sig/DiffBind_results_K27ac_1000_up_peaks_significant.csv")
k27ac_down <- read.csv("~/Agastya_UM/M_ChipSeq/M_OP/batch2_t2p/DiffBind/summits_1000/K27ac/sig/DiffBind_results_K27ac_1000_down_peaks_significant.csv")

k27me3_up   <- read.csv("~/Agastya_UM/M_ChipSeq/M_OP/batch2_t2p/DiffBind/summits_1000/K27me3/sig/DiffBind_results_K27me3_1000_up_peaks_significant.csv")
k27me3_down <- read.csv("~/Agastya_UM/M_ChipSeq/M_OP/batch2_t2p/DiffBind/summits_1000/K27me3/sig/DiffBind_results_K27me3_1000_down_peaks_significant.csv")

# --------------------------
# Combine data
# --------------------------

df <- rbind(
  data.frame(log2FC = k27ac_up$Fold,    Mark = "H3K27ac"),
  data.frame(log2FC = k27ac_down$Fold,  Mark = "H3K27ac"),
  data.frame(log2FC = k27me3_up$Fold,   Mark = "H3K27me3"),
  data.frame(log2FC = k27me3_down$Fold, Mark = "H3K27me3")
)

df <- na.omit(df)

# --------------------------
# Sample sizes
# --------------------------

n_ac  <- sum(df$Mark == "H3K27ac")
n_me3 <- sum(df$Mark == "H3K27me3")

label_text <- paste0(
  "n(H3K27ac) = ", format(n_ac, big.mark = ","),
  "\n",
  "n(H3K27me3) = ", format(n_me3, big.mark = ",")
)

# --------------------------
# Plot
# --------------------------

ymax <- max(df$log2FC, na.rm = TRUE)
ymin <- min(df$log2FC, na.rm = TRUE)

p <- ggplot(df, aes(x = Mark, y = log2FC, fill = Mark)) +
  
  geom_violin(
    trim = FALSE,
    width = 1.6,
    alpha = 0.60,
    color = NA
  ) +
  
  geom_boxplot(
    width = 0.35,
    outlier.shape = NA,
    linewidth = 0.4,
    color = "black",
    alpha = 0.9
  ) +
  
  geom_jitter(
    width = 0.12,
    size = 0.35,
    alpha = 0.12,
    color = "black"
  ) +
  
  scale_x_discrete(
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.5
  ) +
  
  annotate(
    "label",
    x = 2.45,
    y = ymax * 0.95,
    label = label_text,
    hjust = 1,
    vjust = 1,
    size = 4,
    fill = "white",
    label.size = 0.5
  ) +
  
  scale_fill_manual(
    values = c(
      "H3K27ac" = "#D55E00",
      "H3K27me3" = "#0072B2"
    )
  ) +
  
  labs(
    x = NULL,
    y = expression(Log[2]~Fold~Change),
    title = "Distribution of Differential Peak Fold Changes"
  ) +
  
  theme_bw(base_size = 14) +
  
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    
    panel.border = element_rect(
      colour = "black",
      linewidth = 0.8
    ),
    
    axis.line = element_blank(),
    
    axis.title.y = element_text(
      size = 16,
      face = "bold"
    ),
    
    axis.text.x = element_text(
      size = 14,
      face = "bold",
      colour = "black"
    ),
    
    axis.text.y = element_text(
      size = 12,
      colour = "black"
    ),
    
    axis.ticks = element_line(
      linewidth = 0.8
    ),
    
    legend.position = "none",
    
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 16
    )
  ) +
  
  coord_cartesian(
    ylim = c(ymin, ymax * 1.1),
    clip = "off"
  )

p

ggsave(
  "K27ac_vs_K27me3_log2FC_publication.pdf",
  p,
  width = 4.5,
  height = 5,
  dpi = 600
)

#---------------------------------------------------------
# Box plot of H3K27ac on H3K27me3 and H3K27me3 on H3K27ac
#---------------------------------------------------------

oppo <- read.csv("~/Agastya_UM/M_ChipSeq/M_OP/batch2_t2p/DiffBind/summits_1000/K27me3/original_coordinates/oc/5k/oppo_diffbind/me3Peaks_acSignals_DiffBind_results.csv") 
k27ac_up <- oppo[oppo$Fold>0,]
k27ac_down <- oppo[oppo$Fold<0,]
#k27ac_up   <- read.csv("~/Agastya_UM//M_ChipSeq/M_OP/batch2_t2p/DiffBind/summits_1000/K27ac/sig/DiffBind_results_K27ac_1000_up_peaks_significant.csv", header = TRUE)
#k27ac_down <- read.csv("~/Agastya_UM/M_ChipSeq/M_OP/batch2_t2p/DiffBind/summits_1000/K27ac/sig/DiffBind_results_K27ac_1000_down_peaks_significant.csv", header = TRUE)

k27me3_up   <- read.csv("~/Agastya_UM/M_ChipSeq/M_OP/batch2_t2p/DiffBind/summits_1000/K27me3/sig/DiffBind_results_K27me3_1000_up_peaks_significant.csv", header = TRUE)
k27me3_down <- read.csv("~/Agastya_UM/M_ChipSeq/M_OP/batch2_t2p/DiffBind/summits_1000/K27me3/sig/DiffBind_results_K27me3_1000_down_peaks_significant.csv", header = TRUE)

k27ac_up$Mark   <- "H3K27ac"
k27ac_down$Mark <- "H3K27ac"

k27me3_up$Mark   <- "H3K27me3"
k27me3_down$Mark <- "H3K27me3"

df <- rbind(
  data.frame(log2FC = k27ac_up$Fold,   Mark = "H3K27ac"),
  #data.frame(log2FC = k27ac_down$Fold, Mark = "H3K27ac"),
  data.frame(log2FC = k27me3_up$Fold,   Mark = "H3K27me3"),
  data.frame(log2FC = k27me3_down$Fold, Mark = "H3K27me3")
)

library(ggplot2)

p <- ggplot(df, aes(x = Mark, y = log2FC, fill = Mark)) +
  
  geom_jitter(
    width = 0.15,
    size = 0.5,
    alpha = 0.35
  ) +
  
  geom_boxplot(
    outlier.shape = NA,
    width = 0.6,
    alpha = 0.8
  ) +
  
  geom_hline(yintercept = 0, linetype = "dashed") +
  
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  labs(
    title = "Differential peak log2FC distribution across H3K27me3 Peaks",
    y = "log2 Fold Change",
    x = ""
  ) +
  scale_fill_manual(
    values = c("H3K27ac" = "#D55E00", "H3K27me3" = "#0072B2")
  ) +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5)
  )

ggsave(
  "K27me3_K27ac_vs_K27me3_log2FC_boxplot.pdf",
  p,
  width = 6,
  height = 4
)

p


