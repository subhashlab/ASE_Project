library(DiffBind)

dbObj <- dba(sampleSheet = "/iitjfs/home/2025rbb1002/M_ChipSeq/M_OP/Samples_files/SampleK27me3_24.csv")

dbObj <- dba.count(
  dbObj,
  summits = 1000,
  score = DBA_SCORE_READS
)

dba.show(dbObj)

dbObj <- dba.contrast(
  dbObj,
  categories = DBA_CONDITION,
  minMembers = 2
)

dba.show(dbObj, bContrasts = TRUE)

dbObj <- dba.analyze(
  dbObj,
  method = DBA_DESEQ2
)

diff_peaks <- dba.report(dbObj)

diff_peaks

diff_df <- as.data.frame(diff_peaks)
write.csv(diff_df, "DiffBind_results_K27me3_1000.csv", row.names = FALSE)
