module load bioinfo/gcc-13.1.0
module load bioinfo/R-4.5.1
export LD_LIBRARY_PATH=/iitjfs/apps/compiler/gcc-13.1.0/lib64
R --vanilla


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
library(tidyverse)
library(ggrepel)

#---------------------------
# Read annotated peaks
#---------------------------
peak_df <- read.csv("/home/2024rbb1003/projects/ASE_NB_project/eCLIP/results/eclip_peaks2_annotated.csv")

#---------------------------
# Thresholds
#---------------------------
p_cutoff <- -log10(0.05)      # 1.30103
signal_cutoff <- 3

#---------------------------
# Categorize peaks
#---------------------------
peak_df <- peak_df %>%
  mutate(
    category = case_when(
      signalValue >= signal_cutoff &
        pValue >= p_cutoff ~ "Significant",
      TRUE ~ "Not Significant"
    ),

    HDAC7 = ifelse(SYMBOL == "HDAC7",
                   "HDAC7",
                   "Other")
  )

#---------------------------
# HDAC7 peaks for labels
#---------------------------
# Only significant HDAC7 peaks
hdac7_peaks <- peak_df %>%
  filter(
    SYMBOL == "HDAC7",
    signalValue >= signal_cutoff,
    pValue >= p_cutoff
  )

#---------------------------
# Volcano plot
#---------------------------
volcano <- ggplot() +

  ## all non-HDAC7 peaks
  geom_point(
    data = peak_df %>% filter(SYMBOL != "HDAC7"),
    aes(signalValue, pValue, fill = category),
    shape = 21,
    color = "black",
    stroke = 0.25,
    size = 2,
    alpha = 0.7
  ) +

  ## HDAC7 peaks
  geom_point(
    data = hdac7_peaks,
    aes(signalValue, pValue),
    shape = 21,
    fill = "#d73027",
    color = "black",
    stroke = 0.6,
    size = 4
  ) +

  ## threshold lines
  geom_hline(
    yintercept = p_cutoff,
    linetype = "dashed",
    color = "grey40"
  ) +

  geom_vline(
    xintercept = signal_cutoff,
    linetype = "dashed",
    color = "grey40"
  ) +

  ## label HDAC7 peaks
  geom_text_repel(
    data = hdac7_peaks,
    aes(signalValue, pValue,
        label = paste0(annotation)),
    size = 3.5,
    max.overlaps = Inf
  ) +

  scale_fill_manual(values = c(
    "Significant" = "#2C7BB6",
    "Not Significant" = "grey85"
  )) +

  labs(
    title = "HNRNPK eCLIP Peaks",
    x = "Signal Value",
    y = expression(-log[10](italic(P))),
    fill = ""
  ) +

  theme_bw(base_size = 13) +

  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "top"
  )

volcano

ggsave(
  "HNRNPK_eCLIP_volcano_HDAC7.pdf",
  volcano,
  width = 6,
  height = 6
)

volcano <- ggplot() +

  ## all non-HDAC7 peaks
  geom_point(
    data = peak_df %>% filter(SYMBOL != "HDAC7"),
    aes(signalValue, pValue, fill = category),
    shape = 21,
    color = "black",
    stroke = 0.25,
    size = 2,
    alpha = 0.7
  ) +

  ## HDAC7 peaks
  geom_point(
    data = hdac7_peaks,
    aes(signalValue, pValue),
    shape = 21,
    fill = "#d73027",
    color = "black",
    stroke = 0.6,
    size = 4
  ) +

  ## threshold lines
  geom_hline(
    yintercept = p_cutoff,
    linetype = "dashed",
    color = "grey40"
  ) +

  geom_vline(
    xintercept = signal_cutoff,
    linetype = "dashed",
    color = "grey40"
  ) +

  scale_fill_manual(values = c(
    "Significant" = "#2C7BB6",
    "Not Significant" = "grey85"
  )) +

  labs(
    title = "HNRNPK eCLIP Peaks",
    x = "Signal Value",
    y = expression(-log[10](italic(P))),
    fill = ""
  ) +

  theme_bw(base_size = 13) +

  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "top"
  )


ggsave(
  "HNRNPK_eCLIP_volcano_HDAC7_2.pdf",
  volcano,
  width = 6,
  height = 6
)



