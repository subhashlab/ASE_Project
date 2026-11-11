library(DiffBind)

db <- dba(sampleSheet="/iitjfs/home/2025rbb1002/M_ChipSeq/M_OP/batch2_t2p/DiffBind/summits_1000/K27ac/original_coordinates/oc/5k/oppo_diffbind/oppo_acPeaks_me3Signals.csv")
db <- dba.count(db, summits=FALSE, filter=0)   # VERY important
db <- dba.contrast(db, categories=DBA_CONDITION, minMembers=2)
db <- dba.analyze(db)
res <- dba.report(db, th = 1)
res <- dba.report(db)
res <- dba.report(db, th = 1, bCounts = TRUE)

write.csv(as.data.frame(res),
          file = "/iitjfs/home/2025rbb1002/M_ChipSeq/M_OP/batch2_t2p/DiffBind/summits_1000/K27ac/original_coordinates/oc/5k/oppo_diffbind/acPeaks_me3Signals_DiffBind_results.csv",
          row.names = FALSE)
          
write.table(as.data.frame(res),
            file = "/iitjfs/home/2025rbb1002/M_ChipSeq/M_OP/batch2_t2p/DiffBind/summits_1000/K27ac/original_coordinates/oc/5k/oppo_diffbind/acPeaks_me3Signals_DiffBind_results.txt",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)
            
saveRDS(db, file = "/iitjfs/home/2025rbb1002/M_ChipSeq/M_OP/batch2_t2p/DiffBind/summits_1000/K27ac/original_coordinates/oc/5k/oppo_diffbind/acPeaks_me3Signals_DiffBind_results.rds")

df <- data.frame(res)

df$category <- "NotSig"
df$category[df$FDR < 0.05 & df$Fold > 0.5] <- "Enriched"
df$category[df$FDR < 0.05 & df$Fold < -0.5] <- "Depleted"

#Plot
library(ggplot2)
library(ggrepel)
library(dplyr)

p <- ggplot(df, aes(x = Fold, y = -log10(FDR), color = category)) +
  geom_point(size = 2) +
  scale_color_manual(values = c("Depleted" = "blue",
                                "NotSig" = "grey",
                                "Enriched" = "red")) +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  theme_minimal() +
  labs(
    title = "acPeaks_me3Signals Volcano Plot",
    x = "log2 Fold Change",
    y = "-log10(FDR)"
  )

ggsave("/iitjfs/home/2025rbb1002/M_ChipSeq/M_OP/batch2_t2p/DiffBind/summits_1000/K27ac/original_coordinates/oc/5k/oppo_diffbind/acPeaks_me3Signals_DiffBind_volcano.pdf", width=6, height=6)

