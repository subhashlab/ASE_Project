library(dplyr)
library(readr)
 
res_all <- read_csv(
   "/home/2024rbb1003/projects/ASE_NB_project/RIPSeq/transcript_level_results/edgeR_ALL_contrasts_ALL_transcripts_TEs_annotated.csv",
   show_col_types = FALSE
)
                                                                                                                  
log2fc_cutoff <- 1      # ≥ 2-fold enrichment
p_cutoff      <- 0.05
final_candidates <- res_all %>%
  filter(
     sele_2_log2FC >= log2fc_cutoff,           # dsRIP enrichment (direction)
     dsRIP_p_val <= p_cutoff,                  # statistical support in dsRIP
     (
       INPUT_p_val > p_cutoff |                # no significant INPUT change
       abs(INPUT_log2FC) < log2fc_cutoff        # or weak INPUT change
     )
  )
write.csv(
   final_candidates,
   "dsRIP_specific_enriched_candidates.csv",
   row.names = FALSE
)
nrow(res_all)            # total tested genes
#[1] 5079170
nrow(final_candidates)  # final candidates
#[1] 18130
 
# How many fail each filter
table(res_all$sele_2_log2FC >= log2fc_cutoff)

#  FALSE    TRUE 
#4195062  884108 
table(res_all$dsRIP_p_val <= p_cutoff)

#  FALSE    TRUE 
#5022459   56711  
 
transcript_info <- read.csv("/iitjfs/home/2024rbb1003/projects/ASE_NB_project/RIPSeq/TE_annotation.csv", header = TRUE, stringsAsFactors = FALSE)
res_df <- final_candidates %>%
  left_join(transcript_info, by = "transcript_id")
write.csv(res_df, "dsRIP_specific_enriched_candidates_TE_annotated.csv", row.names = FALSE)
