#DiffBind
library(DiffBind)

db <- dba(sampleSheet="/iitjfs/home/2025rbb1002/M_ChipSeq/M_OP/batch2_t2p/DiffBind/summits_1000/K27ac/original_coordinates/oc/len/SE_MChip/diffbind_SE/Union/SE_sample.csv")
db <- dba.count(db, summits=FALSE)   # VERY important
db <- dba.contrast(db, categories=DBA_CONDITION, minMembers=2)
db <- dba.analyze(db)

res <- dba.report(db)

write.csv(as.data.frame(res),
          file = "/iitjfs/home/2025rbb1002/M_ChipSeq/M_OP/batch2_t2p/DiffBind/summits_1000/K27ac/original_coordinates/oc/len/SE_MChip/diffbind_SE/Union/SE_DiffBind_results.csv",
          row.names = FALSE)
          
write.table(as.data.frame(res),
            file = "/iitjfs/home/2025rbb1002/M_ChipSeq/M_OP/batch2_t2p/DiffBind/summits_1000/K27ac/original_coordinates/oc/len/SE_MChip/diffbind_SE/Union/SE_DiffBind_results.txt",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)
            
saveRDS(db, file = "/iitjfs/home/2025rbb1002/M_ChipSeq/M_OP/batch2_t2p/DiffBind/summits_1000/K27ac/original_coordinates/oc/len/SE_MChip/diffbind_SE/Union/diffbind_SE.rds")


res_sig <- dba.report(db, th=0.05)

write.csv(as.data.frame(res_sig),
          file = "/iitjfs/home/2025rbb1002/M_ChipSeq/M_OP/batch2_t2p/DiffBind/summits_1000/K27ac/original_coordinates/oc/len/SE_MChip/diffbind_SE/Union/SE_significant.csv",
          row.names = FALSE)
          
pdf("/iitjfs/home/2025rbb1002/M_ChipSeq/M_OP/batch2_t2p/DiffBind/summits_1000/K27ac/original_coordinates/oc/len/SE_MChip/diffbind_SE/Union/PCA_SE.pdf")
dba.plotPCA(db, DBA_CONDITION)
dev.off()

pdf("/iitjfs/home/2025rbb1002/M_ChipSeq/M_OP/batch2_t2p/DiffBind/summits_1000/K27ac/original_coordinates/oc/len/SE_MChip/diffbind_SE/Union/heatmap_SE.pdf")
dba.plotHeatmap(db)
dev.off()

#Naming the SEs
se <- read.table("/iitjfs/home/2025rbb1002/M_ChipSeq/M_OP/batch2_t2p/DiffBind/summits_1000/K27ac/original_coordinates/oc/len/SE_MChip/diffbind_SE/Union/SE_union_named.bed", header=FALSE)
colnames(se) <- c("Chr","Start","End","SE_name")

df <- read.csv("/iitjfs/home/2025rbb1002/M_ChipSeq/M_OP/batch2_t2p/DiffBind/summits_1000/K27ac/original_coordinates/oc/len/SE_MChip/diffbind_SE/Union/SE_significant.csv")
colnames(df) <- c("Chr","Start","End","width", "strand", "Conc", "Conc_HNRNPKsi", "Conc_Control", "Fold", "p.value", "FDR")

df <- merge(df, se,
            by=c("Chr","Start","End"),
            all.x=TRUE)

df$category <- "NotSig"
df$category[df$FDR < 0.05 & df$Fold > 0.25] <- "Enriched"
df$category[df$FDR < 0.05 & df$Fold < -0.25] <- "Depleted"

#Plot
library(ggplot2)
library(ggrepel)
library(dplyr)


# Top 10 upregulated SEs
top_up <- df %>%
  filter(FDR < 0.05, Fold > 0.25) %>%
  arrange(desc(Fold)) %>%
  slice_head(n = 10)

# Top 10 downregulated SEs
top_down <- df %>%
  filter(FDR < 0.05, Fold < -0.25) %>%
  arrange(Fold) %>%
  slice_head(n = 10)

top_labels <- bind_rows(top_up, top_down)

p <- ggplot(df, aes(x = Fold, y = -log10(FDR), color = category)) +
  geom_point(size = 2) +
  geom_text_repel(
    data = top_labels,
    aes(label = SE_name),   # replace SE_name with your column name
    size = 3,
    max.overlaps = Inf
  ) +
  scale_color_manual(values = c("Depleted" = "blue",
                                "NotSig" = "grey",
                                "Enriched" = "red")) +
  geom_vline(xintercept = c(-0.25, 0.25), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  theme_minimal() +
  labs(
    title = "Super-enhancer Volcano Plot",
    x = "log2 Fold Change",
    y = "-log10(FDR)"
  )

ggsave("/iitjfs/home/2025rbb1002/M_ChipSeq/M_OP/batch2_t2p/DiffBind/summits_1000/K27ac/original_coordinates/oc/len/SE_MChip/diffbind_SE/Union/SE_volcano_ggplot_0.25_named.pdf", width=6, height=6)

library(ggplot2)

df <- as.data.frame(res)

# Categorize
df$category <- "NotSig"
df$category[df$FDR < 0.05 & df$Fold > 0.25] <- "Up"
df$category[df$FDR < 0.05 & df$Fold < -0.25] <- "Down"

# Save results
sig <- df[df$FDR < 0.05 & abs(df$Fold) > 0.25,]
up  <- df[df$FDR < 0.05 & df$Fold > 0.25,]
down<- df[df$FDR < 0.05 & df$Fold < -0.25,]

write.csv(sig,  "SE_sig_FC0.25_FDR0.05.csv", row.names=FALSE)
write.csv(up,   "SE_upregulated.csv", row.names=FALSE)
write.csv(down, "SE_downregulated.csv", row.names=FALSE)

# Volcano plot with cutoff lines
p <- ggplot(df, aes(x=Fold, y=-log10(FDR), color=category)) +
  geom_point(size=2) +
  scale_color_manual(values=c("blue","grey","red")) +
  
  # Fold-change cutoffs
  geom_vline(xintercept=c(-0.25,0.25), linetype="dashed", color="black") +
  
  # P-value cutoff (THIS is your abline equivalent)
  geom_hline(yintercept=-log10(0.05), linetype="dashed", color="black") +
  
  theme_minimal() +
  labs(title="Super-enhancer Volcano Plot",
       x="log2 Fold Change",
       y="-log10(FDR)")

ggsave("/iitjfs/home/2025rbb1002/M_ChipSeq/M_OP/batch2_t2p/DiffBind/summits_1000/K27ac/original_coordinates/oc/len/SE_MChip/diffbind_SE/Union/SE_volcano_FC0.25_FDR0.05.pdf", plot=p, width=6, height=6)

# Count Up and Down
up_count   <- sum(df$FDR < 0.05 & df$Fold > 0.25)
down_count <- sum(df$FDR < 0.05 & df$Fold < -0.25)

up_count = 36
down_count = 19

up <- df[df$FDR < 0.05 & df$Fold > 0.25, ]
down <- df[df$FDR < 0.05 & df$Fold < -0.25, ]

write.table(
  up[, c("Chr","Start","End","SE_name")],
  "/iitjfs/home/2025rbb1002/M_ChipSeq/M_OP/batch2_t2p/DiffBind/summits_1000/K27ac/original_coordinates/oc/len/SE_MChip/diffbind_SE/Union/SE_upregulated.bed",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

write.table(
  down[, c("Chr","Start","End","SE_name")],
  "/iitjfs/home/2025rbb1002/M_ChipSeq/M_OP/batch2_t2p/DiffBind/summits_1000/K27ac/original_coordinates/oc/len/SE_MChip/diffbind_SE/Union/SE_downregulated.bed",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)
