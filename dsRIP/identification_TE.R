library(limma)
library(edgeR)
library(tidyverse)
library(readxl)

# Load filtered counts
counts_df <- read.delim("/iitjfs/home/2024rbb1003/projects/ASE_NB_project/RIPSeq/quantification/hg38_counts_TE_cleaned.txt", row.names = 1, check.names = FALSE)

# Load metadata
meta <- read_excel("/iitjfs/home/2024rbb1003/projects/ASE_NB_project/RIPSeq/metadata_RIP.xlsx")

# Ensure Sample_ID in metadata matches column names in counts
counts_df <- counts_df[, meta$Sample_ID]

# Create a combined group variable for edgeR
meta <- meta %>%
  mutate(group = paste(fraction, condition, sep = "_"))

# Create DGEList
dge <- DGEList(counts = counts_df, group = meta$group, samples = meta)

# Calculate normalization factors using TMM
dge <- calcNormFactors(dge, method = "TMM")

# Optional: compute normalized library sizes
norm_factors <- (dge$samples$lib.size * dge$samples$norm.factors) / mean(dge$samples$lib.size * dge$samples$norm.factors)
names(norm_factors) <- dge$samples$Sample_ID

# Save normalization factors for reference
norm_df <- data.frame(sample = names(norm_factors), div_norm = norm_factors)
write.csv(norm_df, file = "dge_normalization_factors.csv", row.names = FALSE)



# log2 CPM values (normalized)
logCPM_mat <- cpm(dge, log = TRUE)

# Optional: add gene annotation if you have it
# r_anno should contain transcript/gene info
# rownames(logCPM_mat) must match rownames(r_anno)

logCPM_df <- as.data.frame(logCPM_mat) %>%
  mutate(id = rownames(logCPM_mat)) %>%
  tidyr::pivot_longer(-id, names_to = "Sample_ID", values_to = "logCPM") %>%
  left_join(meta, by = "Sample_ID")

library(ggplot2)
library(ggrepel)

# Transpose for sample distance calculation
sample_scaled <- scale(t(logCPM_mat))
sample_dist <- dist(sample_scaled)

# Classical MDS
mds_df <- as.data.frame(cmdscale(sample_dist))
colnames(mds_df) <- c("MDS1", "MDS2")
mds_df$Sample_ID <- rownames(mds_df)

# Merge with metadata for coloring
mds_df <- left_join(mds_df, meta, by = "Sample_ID")

# Plot
mds_plot <- ggplot(mds_df, aes(MDS1, MDS2, color = group)) +
  geom_point(size = 3) +
  geom_text_repel(aes(label = Sample_ID), size = 3) +
  theme_bw() +
  labs(title = "MDS plot of dsRIP samples") +
  theme(legend.position = "right")

pdf("dsRIP_MDS_plot.pdf", width = 5, height = 4)
print(mds_plot)
dev.off()

### Step 34: design, dispersion, GLM fit (quasi-likelihood)
library(edgeR)

# 1) ensure group column exists in dge$samples (it should if you created dge with samples=meta)
# If your samples df column name is 'group' use that; if it's named differently (e.g., 'group'), adjust below.
if(!"group" %in% colnames(dge$samples)){
  stop("dge$samples does not contain a 'group' column. Make sure you created dge with samples=meta and meta has 'group'.")
}

# 2) Build design matrix with one column per group (no intercept) - matches your original script
design <- model.matrix(~0 + group, data = dge$samples)
colnames(design) <- levels(factor(dge$samples$group))
rownames(design) <- rownames(dge$samples)

# Save design to file for record
write.table(design, file = "design_step34.txt", sep = "\t", quote = FALSE, col.names = NA)

# 3) Estimate dispersions using the design (robust = TRUE as in your script)
dge_disp <- estimateDisp(dge, design = design, robust = TRUE)

# Save dispersion summary
disp_summary <- data.frame(
  sample = rownames(dge_disp$samples),
  lib.size = dge_disp$samples$lib.size,
  norm.factors = dge_disp$samples$norm.factors
)
write.csv(disp_summary, file = "dispersion_sample_summary_step34.csv", row.names = FALSE)

# Optional: produce BCV plot and save as PDF (safe device)
pdf("BCV_plot_step34.pdf", width = 5, height = 4)
plotBCV(dge_disp)
dev.off()

# 4) Fit quasi-likelihood GLM (glmQLFit) with robust = TRUE as in your script
dge_fit <- glmQLFit(dge_disp, design = design, robust = TRUE)

# Save fitted object and intermediate objects
save(dge, dge_disp, dge_fit, design, file = "edgeR_objects_step34.RData")

# Quick message
cat("Step 34 complete: design, dispersion estimated, glmQLFit done.\n")
cat("Saved: design_step34.txt, dispersion_sample_summary_step34.csv, BCV_plot_step34.pdf, edgeR_objects_step34.RData\n")


print(colnames(design))
#[1] "dsRIP_ctrl" "dsRIP_trt"  "INPUT_ctrl" "INPUT_trt" 

#r_anno <- read_tsv("/iitjfs/home/2024rbb1003/projects/projectRseq/EndoOrganoids/scripts_combined/gene_info_hg38_113.tsv", show_col_types = FALSE)



# ===============================
# 1️⃣ Assign edgeR matrices
# ===============================
edge_f <- dge$counts           # filtered counts matrix
edge_n <- logCPM_mat           # log2 CPM normalized values

# ===============================
# 2️⃣ Design matrix (already created)
# ===============================
# design <- design  # already matches reference script

# ===============================
# 3️⃣ Gene annotation
# ===============================
# Already loaded:
# r_anno <- read_tsv("gene_info_hg38_113.tsv")

# ===============================
# 4️⃣ Plotting colors
# ===============================
dsRIP_ctrl_col <- "#1f77b4"
dsRIP_trt_col  <- "#ff7f0e"
input_trt_col   <- "#2ca02c"

# ===============================
# 5️⃣ Quick sanity check
# ===============================
stopifnot(exists("edge_f"), exists("edge_n"), exists("design"), exists("r_anno"))




edgeRglmQLF<-function(mat=edge_f, # object of class DGEGLM
                       contro, # comparison, create with makeContrasts
                      cpm_mat=edge_n, #used to calculate average signal
                       label="", # include in label column
                       sig_thr=0.5, # signal threshold
                       sig_col="CPM", # signal column
                       fc_thr=0.5, # fold change threshold
                       pval_thr=0.05, # p-value threshold
                       pval_col="p_val", # p-value column
                       names=FALSE){ # include label in column names
   degs<-glmQLFTest(edge_f,contrast=contro)$table[,-3]
   # degs<-glmTreat(edge_f,contrast=contro,lfc=log2(1.2))$table[,-3] # alternative with FC threshold
   colnames(degs)<-c("log2_FC","log2_CPM","p_val")
   a_levels<-rownames(contro)[which(contro!=0)]
   a_samples<-which(cpm_mat$samples$group%in%a_levels)
   cpm_sele<-cpm(cpm_mat,log=T)[,a_samples]
  degs$log2_CPM<-apply(cpm_sele,1,function(x) mean(x))
   #degs<-exactTest(edge_c, pair=cond, dispersion=bcv^2)$table
   degs$p_adj<-p.adjust(degs$p_val, method ="BH")
   degs$CPM<-round(2^degs$log2_CPM,2)
   degs$class<-"="
   degs[which(degs[,sig_col]>=sig_thr & degs$log2_FC>=fc_thr & degs[,pval_col]<=pval_thr),"class"]<-"+"
   degs[which(degs[,sig_col]>=sig_thr & degs$log2_FC<=(-fc_thr) & degs[,pval_col]<=pval_thr),"class"]<-"-"
   degs$class<-as.factor(degs$class)
   degs$comp<-label
   degs$id<-rownames(degs)
   degs<-degs[,c("id","comp","log2_FC","CPM","p_val","p_adj","class")]
   if(names=="TRUE"){
     newnames<-paste(label,colnames(degs),sep="_")
     colnames(degs)<-newnames
   }
   return(degs)
 }
 
contro_label <- "(dsRIP_trt-dsRIP_ctrl)"
contro <- makeContrasts(contro_label,levels=design) 

> head(counts_df)
                BK846.002L0001 BK846.002L0002 BK846.002L0003 BK846.002L0004
ENSG00000157911            360             58            175            239
ENSG00000142655            680            246            469            563
ENSG00000149527            155             54             48             82
ENSG00000171621            602            325            358            473
ENSG00000173614            353            100             43            271
ENSG00000215720             21             36             21             20
                BK846.002L0005 BK846.002L0006 BK846.01L0013 BK846.01L0014
ENSG00000157911            231            234           315            56
ENSG00000142655            674            698           308           170
ENSG00000149527            215            110            18             6
ENSG00000171621            416            508           619           233
ENSG00000173614            150             35           867           138
ENSG00000215720             36             29             2             2
                BK846.01L0015 BK846.01L0016 BK846.01L0017 BK846.01L0018
ENSG00000157911           144            96           171           112
ENSG00000142655           572            93           454           401
ENSG00000149527            18             0            24            14
ENSG00000171621           728           140           399           497
ENSG00000173614           227           569           562           192
ENSG00000215720             2             0             0             0
> 
> 
> 
> 

> head(meta)
# A tibble: 6 × 14
  sample_names         Sample_ID fraction condition Condition treatment cellline
  <chr>                <chr>     <chr>    <chr>     <chr>     <chr>     <chr>   
1 Unknown_BK846-002L0… BK846.00… INPUT    ctrl      Control   Csi       HeLa    
2 Unknown_BK846-002L0… BK846.00… INPUT    trt       HNRNPKsi  siRNA1    HeLa    
3 Unknown_BK846-002L0… BK846.00… INPUT    trt       HNRNPKsi  siRNA1    HeLa    
4 Unknown_BK846-002L0… BK846.00… INPUT    ctrl      Control   Csi       HeLa    
5 Unknown_BK846-002L0… BK846.00… INPUT    trt       HNRNPKsi  siRNA2    HeLa    
6 Unknown_BK846-002L0… BK846.00… INPUT    trt       HNRNPKsi  siRNA2    HeLa    
# ℹ 7 more variables: other_information <chr>, extraction <chr>, batch <chr>,
#   library <chr>, BAM_file_location <chr>, FASTQ_file_location <chr>,
#   group <chr>


# Assign variables expected by edgeRglmQLF
edge_f <- dge_fit
edge_n <- dge
contro_label <- "(dsRIP_trt-dsRIP_ctrl)"
contro <- makeContrasts(contro_label,levels=design) 
dsRIP_DEGs<-edgeRglmQLF(mat=edge_f, cpm_mat=edge_n, contro=contro, label="dsRIP", sig_thr=0.5, sig_col="CPM", fc_thr=0, pval_thr=0.05, pval_col="p_val",names=T)
 
contro_label <- "(INPUT_trt-INPUT_ctrl)"
contro <- makeContrasts(contro_label,levels=design) 
INPUT_DEGs<-edgeRglmQLF(mat=edge_f, cpm_mat=edge_n, contro=contro, label="INPUT", sig_thr=0.5, sig_col="CPM", fc_thr=0, pval_thr=0.05, pval_col="p_val",names=T)
 
contro_label <- "(dsRIP_trt + dsRIP_ctrl) - (INPUT_trt + INPUT_ctrl)"
contro <- makeContrasts(contro_label,levels=design) 
sele1_DEGs<-edgeRglmQLF(mat=edge_f, cpm_mat=edge_n, contro=contro, label="sele_1", sig_thr=0.5, sig_col="CPM", fc_thr=0, pval_thr=0.05, pval_col="p_val",names=T)
 
### Selection 2
contro_label <- "(dsRIP_trt - dsRIP_ctrl) - (INPUT_trt - INPUT_ctrl)"
contro <- makeContrasts(contro_label,levels=design) 
sele2_DEGs<-edgeRglmQLF(mat=edge_f, cpm_mat=edge_n, contro=contro, label="sele_2", sig_thr=0.5, sig_col="CPM", fc_thr=0, pval_thr=0.05, pval_col="p_val",names=T)
dsRIP_1f_df<- bind_cols(sele1_DEGs[,c(1,3,4,5,7)], sele2_DEGs[,c(3,5,7)], dsRIP_DEGs[,c(3,5,7)], INPUT_DEGs[,c(3,5,7)])
 

#fuse_df<-full_join(r_anno , dsRIP_1f_df, by=c("gene_name" ="sele_1_id"))
fuse_df<- dsRIP_1f_df
fuse_df$sele_3_class <- ifelse(abs(fuse_df$dsRIP_log2_FC)>abs(fuse_df$INPUT_log2_FC),"+","=")
 
fuse_df$dsRIP_final_class<-"=" # Final identification of J2 candidates
 
fuse_df[which(fuse_df$sele_1_class=="+" &
              fuse_df$sele_2_class=="+" &
              fuse_df$sele_3_class=="+"
              ),"dsRIP_final_class"]<-"+"


fuse_df$dsRIP_final_class <- factor(fuse_df$dsRIP_final_class,levels=c("=","+"))

colnames(fuse_df)
 [1] "gene_id"           "gene_name"         "chromosome"       
 [4] "start"             "end"               "gene_biotype"     
 [7] "sele_1_log2_FC"    "sele_1_CPM"        "sele_1_p_val"     
[10] "sele_1_class"      "sele_2_log2_FC"    "sele_2_p_val"     
[13] "sele_2_class"      "dsRIP_log2_FC"     "dsRIP_p_val"      
[16] "dsRIP_class"       "INPUT_log2_FC"     "INPUT_p_val"      
[19] "INPUT_class"       "sele_3_class"      "dsRIP_final_class"

fuse_df_cands <- fuse_df %>% 
  dplyr::filter(dsRIP_final_class == "+") %>% 
  dplyr::select(dsRIP_log2_FC, INPUT_log2_FC, 
                sele_1_class, sele_2_class, sele_3_class) %>% 
  dplyr::arrange(desc(dsRIP_log2_FC))
in_mat<- fuse_df

in_mat$J2_final_class <- as.factor(in_mat$dsRIP_final_class)

sc<-ggplot(in_mat, aes(x=sele_2_log2_FC,y=sele_1_log2_FC,colour=dsRIP_final_class,alpha=dsRIP_final_class,size=dsRIP_final_class))+
   geom_point(shape=16)+
   geom_hline(yintercept=0,col="grey40",lty=2,lwd=0.6) +
   geom_vline(xintercept=0,col="grey40",lty=2,lwd=0.6) +
   theme_bw() +
   scale_colour_manual(values=c(dsRIP_ctrl_col,dsRIP_trt_col)) +
   scale_size_manual(values=c(0.5,1.5)) +
   scale_alpha_manual(values=c(0.1,0.8))+
   theme(legend.justification=c(0,1), legend.position=c(0,1))+
   theme(legend.title=element_blank())+
   theme(legend.position = "none")+
   theme(panel.grid.minor = element_blank())+
   labs(x = "dsRIP trt vs ctrl Fold Enrichment \n(log2, normalized for INPUT)", y = "dsRIP/INPUT Fold Enrichment \n(log2)")+
   coord_cartesian(xlim = c(-1,3), ylim = c(-3,10))

sc<-sc +
       geom_text_repel(data = subset(in_mat,
                                     dsRIP_final_class=="+" & 
                                     (sele_2_log2_FC^2 + (sele_1_log2_FC/4)^2)>3  # labeling outside a circle
                                     ),
                       aes(x=sele_2_log2_FC,y=sele_1_log2_FC,label=external_gene_name),
                       xlim  = c(1,NA),
                       ylim  = c(2,NA),
                       colour=input_trt_col,
                       alpha=1,
                       fontface=4,
                       size=3,
                       box.padding = unit(0.1, "lines"),
                       point.padding = unit(0.2, "lines"),
                       segment.color = input_trt_col,
                       segment.size=0.1,
                       segment.alpha=0.5) 
 
sc<-ggplot(in_mat, aes(x=sele_2_log2_FC,y=sele_1_log2_FC,colour=dsRIP_final_class,alpha=dsRIP_final_class,size=dsRIP_final_class))+
   geom_point(shape=16)+
   geom_hline(yintercept=0,col="grey40",lty=2,lwd=0.6) +
   geom_vline(xintercept=0,col="grey40",lty=2,lwd=0.6) +
   theme_bw() +
   scale_colour_manual(values=c(dsRIP_ctrl_col,dsRIP_trt_col)) +
   scale_size_manual(values=c(0.5,1.5)) +
   scale_alpha_manual(values=c(0.1,0.8))+
   theme(legend.justification=c(0,1), legend.position=c(0,1))+
   theme(legend.title=element_blank())+
   theme(legend.position = "none")+
   theme(panel.grid.minor = element_blank())+
   labs(x = "dsRIP trt vs ctrl Fold Enrichment \n(log2, normalized for INPUT)", y = "dsRIP/INPUT Fold Enrichment \n(log2)")+
   coord_cartesian(xlim = c(-1,3), ylim = c(-3,10))
 
sc<-sc +
       geom_text_repel(data = subset(in_mat,
                                     dsRIP_final_class=="+" & 
                                     (sele_2_log2_FC^2 + (sele_1_log2_FC/4)^2)>3  # labeling outside a circle
                                     ),
                       aes(x=sele_2_log2_FC,y=sele_1_log2_FC,label=sele_1_id),
                       xlim  = c(1,NA),
                       ylim  = c(2,NA),
                       colour=input_trt_col,
                       alpha=1,
                       fontface=4,
                       size=3,
                       box.padding = unit(0.1, "lines"),
                       point.padding = unit(0.2, "lines"),
                       segment.color = input_trt_col,
                       segment.size=0.1,
                       segment.alpha=0.5) 
 
#ggsave("J2_RNAseq_scatter.pdf",sc,device=cairo_pdf,width = 3, height = 3, units = c("in"))

pdf("J2_RNAseq_scatter.pdf", width = 10, height = 10)
print(sc)
dev.off()
 
# Determine axis limits from your data
xlims <- range(in_mat$sele_2_log2_FC, na.rm = TRUE)
ylims <- range(in_mat$sele_1_log2_FC, na.rm = TRUE)
 

# Add labels for strong candidates (keep original threshold)
sc <- sc +
   geom_text_repel(data = subset(in_mat,
                                 dsRIP_final_class == "+" &
                                 (sele_2_log2_FC^2 + (sele_1_log2_FC/4)^2) > 3),
                   aes(x = sele_2_log2_FC, y = sele_1_log2_FC, label = sele_1_id),
                   xlim  = c(1, NA),
                   ylim  = c(2, NA),
                   colour = input_trt_col,
                   alpha = 1,
                   fontface = 4,
                   size = 3,
                   box.padding = unit(0.1, "lines"),
                   point.padding = unit(0.2, "lines"),
                   segment.color = input_trt_col,
                   segment.size = 0.1,
                   segment.alpha = 0.5)

#ggsave("J2_RNAseq_scatter_2.pdf",sc,device=cairo_pdf,width = 3, height = 3, units = c("in"))
pdf("J2_RNAseq_scatter_2.pdf", width = 10, height = 10)
print(sc)
dev.off()
 
table(is.na(in_mat$sele_1_log2_FC))

#FALSE  TRUE 
#18096 75190 
table(is.na(in_mat$sele_2_log2_FC))

#FALSE  TRUE 
#18096 75190 
table(in_mat$dsRIP_final_class)

#    =     + 
#93279     7 
 
library(ggplot2)
library(ggrepel)
library(dplyr)
 
# Make a copy of the data for plotting
plot_df <- in_mat
 
# Convert dsRIP_final_class to factor for coloring
plot_df$dsRIP_final_class <- factor(plot_df$dsRIP_final_class, levels = c("=", "+"))
 
# Define colors: grey for non-candidates, highlight for positives
cols <- c("=" = "grey70", "+" = "red")  # change "red" to your preferred color
sizes <- c("=" = 0.5, "+" = 2)
alphas <- c("=" = 0.3, "+" = 1)
 
# Scatter plot
sc <- ggplot(plot_df, aes(x = sele_2_log2_FC, y = sele_1_log2_FC)) +
   geom_point(aes(color = dsRIP_final_class, size = dsRIP_final_class, alpha = dsRIP_final_class)) +
   scale_color_manual(values = cols) +
   scale_size_manual(values = sizes) +
   scale_alpha_manual(values = alphas) +
   geom_hline(yintercept = 0, color = "grey40", linetype = 2) +
   geom_vline(xintercept = 0, color = "grey40", linetype = 2) +
   theme_bw() +
   theme(panel.grid.minor = element_blank(),
         legend.title = element_blank(),
         legend.position = "right") +
   labs(x = "dsRIP/INPUT Fold Change (log2, treatment-control)", 
        y = "dsRIP + INPUT Enrichment (log2)") +
   coord_cartesian(xlim = c(min(plot_df$sele_2_log2_FC, na.rm = TRUE)-1, 
                            max(plot_df$sele_2_log2_FC, na.rm = TRUE)+1),
                   ylim = c(min(plot_df$sele_1_log2_FC, na.rm = TRUE)-1, 
                            max(plot_df$sele_1_log2_FC, na.rm = TRUE)+1))
 
# Add labels for positive candidates only
sc <- sc + geom_text_repel(data = subset(plot_df, dsRIP_final_class == "+"),
                            aes(label = sele_1_id),
                            size = 3,
                            fontface = "bold",
                            box.padding = unit(0.2, "lines"),
                            point.padding = unit(0.2, "lines"),
                            segment.size = 0.2,
                            segment.alpha = 0.5)
 

# Save plot
#ggsave("J2_RNAseq_scatter_highlighted.pdf", sc, width = 5, height = 4)
pdf("J2_RNAseq_scatter_highlighted.pdf", width = 10, height = 10)
print(sc)
dev.off()
 
# Save separate list of positive genes
pos_genes <- plot_df %>% filter(dsRIP_final_class == "+") %>%
   select(sele_1_id, dsRIP_log2_FC, INPUT_log2_FC, sele_1_class, sele_2_class, sele_3_class)

write.csv(pos_genes, "J2_RNAseq_candidates.csv", row.names = FALSE)
 
 
 
sc <- sc + geom_text_repel(
   data = subset(plot_df, dsRIP_final_class == "+"),
   aes(label = FALSE),
   size = 3,
   fontface = "bold",
   box.padding = unit(0.2, "lines"),
   point.padding = unit(0.2, "lines"),
   segment.size = 0.2,
   segment.alpha = 0.5,
   max.overlaps = 20   # increase from default
)
#ggsave("J2_RNAseq_scatter_highlighted_2.pdf", sc, width = 5, height = 4)
pdf("J2_RNAseq_scatter_highlighted_2.pdf", width = 10, height = 10)
print(sc)
dev.off()


