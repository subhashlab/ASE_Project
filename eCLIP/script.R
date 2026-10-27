module load bioinfo/gcc-13.1.0
module load bioinfo/R-4.5.1
export LD_LIBRARY_PATH=/iitjfs/apps/compiler/gcc-13.1.0/lib64

zcat GSM2423241_ENCFF241AOZ_peaks_GRCh38.bed.gz \
| awk '$1=="chr12" && $2<47819903 && $3>47782722' \
> HDAC7_rep1.bed

wc -l HDAC7_201_rep1.bed

zcat GSM2423242_ENCFF824IDO_peaks_GRCh38.bed.gz \
| awk '$1=="chr12" && $2<47819903 && $3>47782722' \
> HDAC7_rep2.bed

wc -l HDAC7_201_rep2.bed


zcat GSM2423241_ENCFF241AOZ_peaks_GRCh38.bed.gz \
| awk '$1=="chr12" && $2<47833132 && $3>47782722' \
> HDAC7_rep1.bed

wc -l HDAC7_rep1.bed

zcat GSM2423242_ENCFF824IDO_peaks_GRCh38.bed.gz \
| awk '$1=="chr12" && $2<47833132 && $3>47782722' \
> HDAC7_rep2.bed

wc -l HDAC7_rep2.bed


########################################################################

#input_file <- "eclip_peaks1_annotated.csv"
input_file <- "eclip_peaks2_annotated.csv"

peaks <- read.csv(input_file, stringsAsFactors = FALSE)

cat("Total peaks:", nrow(peaks), "\n")


# Filter HDAC7 intronic peaks

hdac7_introns <- peaks %>%
  filter(
    SYMBOL == "HDAC7",
    grepl("^Intron", annotation)
  )

cat("HDAC7 intronic peaks:", nrow(hdac7_introns), "\n")


# Save all HDAC7 intronic peaks

#write.csv(hdac7_introns, "HDAC7_intronic_peaks_rep1.csv", row.names = FALSE)
write.csv(hdac7_introns, "HDAC7_intronic_peaks_rep2.csv", row.names = FALSE)


# Apply p-value threshold
#
# pValue column contains -log10(p-value)
# p <= 0.05  -->  -log10(p) >= 1.30103

threshold <- -log10(0.05)

hdac7_introns_sig <- hdac7_introns %>%
  filter(pValue >= threshold)

cat("HDAC7 intronic peaks (p <= 0.05):",
    nrow(hdac7_introns_sig), "\n")

#write.csv(hdac7_introns_sig, "HDAC7_intronic_peaks_p0.05_rep1.csv", row.names = FALSE)
write.csv(hdac7_introns_sig, "HDAC7_intronic_peaks_p0.05_rep2.csv", row.names = FALSE)

################################################################################################
library(dplyr)

rep1 <- read.csv("HDAC7_intronic_peaks_p0.05_rep1.csv")
rep2 <- read.csv("HDAC7_intronic_peaks_p0.05_rep2.csv")

# Keep only genomic coordinates
rep1_coords <- rep1 %>%
  select(seqnames, start, end) %>%
  distinct()

rep2_coords <- rep2 %>%
  select(seqnames, start, end) %>%
  distinct()

# Common peaks
common_peaks <- inner_join(
  rep1_coords,
  rep2_coords,
  by = c("seqnames", "start", "end")
)

# Rep1-specific
rep1_only <- anti_join(
  rep1_coords,
  rep2_coords,
  by = c("seqnames", "start", "end")
)

# Rep2-specific
rep2_only <- anti_join(
  rep2_coords,
  rep1_coords,
  by = c("seqnames", "start", "end")
)

cat("Rep1 peaks:", nrow(rep1_coords), "\n")
cat("Rep2 peaks:", nrow(rep2_coords), "\n")
cat("Common:", nrow(common_peaks), "\n")
cat("Rep1 only:", nrow(rep1_only), "\n")
cat("Rep2 only:", nrow(rep2_only), "\n")

write.csv(common_peaks, "Common_HDAC7_peaks.csv", row.names = FALSE)
write.csv(rep1_only, "Rep1_unique_peaks.csv", row.names = FALSE)
write.csv(rep2_only, "Rep2_unique_peaks.csv", row.names = FALSE)

###########################################################################################################################
library(GenomicRanges)
library(dplyr)


rep1 <- read.csv("HDAC7_intronic_peaks_p0.05_rep1.csv",
                 stringsAsFactors = FALSE)

rep2 <- read.csv("HDAC7_intronic_peaks_p0.05_rep2.csv",
                 stringsAsFactors = FALSE)

gr1 <- GRanges(
    seqnames = rep1$seqnames,
    ranges = IRanges(rep1$start, rep1$end),
    strand = rep1$strand
)

gr2 <- GRanges(
    seqnames = rep2$seqnames,
    ranges = IRanges(rep2$start, rep2$end),
    strand = rep2$strand
)

hits <- findOverlaps(gr1, gr2)

cat("Number of overlapping peak pairs:", length(hits), "\n")

idx1 <- queryHits(hits)
idx2 <- subjectHits(hits)

consensus <- data.frame(
    seqnames = rep1$seqnames[idx1],
    start = pmin(rep1$start[idx1], rep2$start[idx2]),
    end   = pmax(rep1$end[idx1], rep2$end[idx2]),
    strand = rep1$strand[idx1],
    SYMBOL = rep1$SYMBOL[idx1],
    annotation = rep1$annotation[idx1],
    rep1_log2FC = rep1$signalValue[idx1],
    rep2_log2FC = rep2$signalValue[idx2],
    rep1_log10P = rep1$pValue[idx1],
    rep2_log10P = rep2$pValue[idx2]
)

cat("Consensus peaks:", nrow(consensus), "\n")

write.csv(
    consensus,
    "HDAC7_consensus_peaks.csv",
    row.names = FALSE
)

bed <- consensus %>%
    transmute(
        chrom = seqnames,
        chromStart = start,
        chromEnd = end,
        name = SYMBOL,
        score = 1000,
        strand = strand
    )

write.table(
    bed,
    "HDAC7_consensus_peaks.bed",
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE
)

cat("Done!\n")


####################################################################################################################
# Create BED6
#peaks <- read.csv("HDAC7_intronic_peaks_rep1.csv", stringsAsFactors = FALSE)
peaks <- read.csv("HDAC7_intronic_peaks_rep2.csv", stringsAsFactors = FALSE)
bed <- peaks %>% transmute( chrom = chrom, chromStart = chromStart, chromEnd = chromEnd, name = name, score = score, strand = strand.1)
#write.table(bed, file = "HDAC7_intronic_peaks_rep1.bed", sep = "\t", quote = FALSE,  row.names = FALSE, col.names = FALSE)
write.table(bed, file = "HDAC7_intronic_peaks_rep2.bed", sep = "\t", quote = FALSE,  row.names = FALSE, col.names = FALSE)



#######################################################################################################################

#input_file <- "eclip_peaks1_annotated.csv"
input_file <- "eclip_peaks2_annotated.csv"

peaks <- read.csv(input_file, stringsAsFactors = FALSE)

cat("Total peaks:", nrow(peaks), "\n")


# Filter HDAC7 intronic peaks

hdac7_introns <- peaks %>%
  filter(
    SYMBOL == "HDAC7")

cat("HDAC7 peaks:", nrow(hdac7_introns), "\n")


# Save all HDAC7 intronic peaks

#write.csv(hdac7_introns, "HDAC7_peaks_rep1.csv", row.names = FALSE)
write.csv(hdac7_introns, "HDAC7_peaks_rep2.csv", row.names = FALSE)


# Apply p-value threshold
#
# pValue column contains -log10(p-value)
# p <= 0.05  -->  -log10(p) >= 1.30103

threshold <- -log10(0.05)

hdac7_introns_sig <- hdac7_introns %>%
  filter(pValue >= threshold)

cat("HDAC7 intronic peaks (p <= 0.05):",
    nrow(hdac7_introns_sig), "\n")

#write.csv(hdac7_introns_sig, "HDAC7_peaks_p0.05_rep1.csv", row.names = FALSE)
write.csv(hdac7_introns_sig, "HDAC7_peaks_p0.05_rep2.csv", row.names = FALSE)

#peaks <- read.csv("HDAC7_peaks_p0.05_rep1.csv", stringsAsFactors = FALSE)
peaks <- read.csv("HDAC7_peaks_p0.05_rep2.csv", stringsAsFactors = FALSE)
bed <- peaks %>% transmute( chrom = chrom, chromStart = chromStart, chromEnd = chromEnd, name = name, score = score, strand = strand.1)
#write.table(bed, file = "HDAC7_peaks_p0.05_rep1.bed", sep = "\t", quote = FALSE,  row.names = FALSE, col.names = FALSE)
write.table(bed, file = "HDAC7_peaks_p0.05_rep2.bed", sep = "\t", quote = FALSE,  row.names = FALSE, col.names = FALSE)




#############################################################################################################################################################################
zcat GSM2423241_ENCFF241AOZ_peaks_GRCh38.bed.gz \
| awk '$1=="chr1" && $2<1268965 && $3>1268696' \
> AluJr_dup74_rep1.bed

wc -l AluJr_dup74_rep1.bed

zcat GSM2423242_ENCFF824IDO_peaks_GRCh38.bed.gz \
| awk '$1=="chr1" && $2<1268965 && $3>1268696' \
> AluJr_dup74_rep2.bed

wc -l AluJr_dup74_rep2.bed



zcat GSM2423241_ENCFF241AOZ_peaks_GRCh38.bed.gz \
| awk '$1=="chr1" && $2<1487504 && $3>1487210' \
> AluY_dup85_rep1.bed

wc -l AluY_dup85_rep1.bed

zcat GSM2423242_ENCFF824IDO_peaks_GRCh38.bed.gz \
| awk '$1=="chr1" && $2<1487504 && $3>1487210' \
> AluY_dup85_rep2.bed

wc -l AluY_dup85_rep2.bed



zcat GSM2423241_ENCFF241AOZ_peaks_GRCh38.bed.gz \
| awk '$1=="chr5" && $2<177488907 && $3>177488606' \
> AluSx1_dup34272_rep1.bed

wc -l AluSx1_dup34272_rep1.bed

zcat GSM2423242_ENCFF824IDO_peaks_GRCh38.bed.gz \
| awk '$1=="chr5" && $2<177488907 && $3>177488606' \
> AluSx1_dup34272_rep2.bed

wc -l AluSx1_dup34272_rep2.bed



zcat GSM2423241_ENCFF241AOZ_peaks_GRCh38.bed.gz \
| awk '$1=="chr5" && $2<177486005 && $3>177485704' \
> AluJr_dup24947_rep1.bed

wc -l AluJr_dup24947_rep1.bed

zcat GSM2423242_ENCFF824IDO_peaks_GRCh38.bed.gz \
| awk '$1=="chr5" && $2<177486005 && $3>177485704' \
> AluJr_dup24947_rep2.bed

wc -l AluJr_dup24947_rep2.bed



zcat GSM2423241_ENCFF241AOZ_peaks_GRCh38.bed.gz \
| awk '$1=="chr5" && $2<177465872 && $3>177465574' \
> AluSx_dup34929_rep1.bed

wc -l AluSx_dup34929_rep1.bed

zcat GSM2423242_ENCFF824IDO_peaks_GRCh38.bed.gz \
| awk '$1=="chr5" && $2<177465872 && $3>177465574' \
> AluSx_dup34929_rep2.bed

wc -l AluSx_dup34929_rep2.bed



###################+400bp
zcat GSM2423241_ENCFF241AOZ_peaks_GRCh38.bed.gz \
| awk '$1=="chr1" && $2<1269365 && $3>1268296' \
> AluJr_dup74_400bp_rep1.bed

zcat GSM2423242_ENCFF824IDO_peaks_GRCh38.bed.gz \
| awk '$1=="chr1" && $2<1269365 && $3>1268296' \
> AluJr_dup74_400bp_rep2.bed

zcat GSM2423241_ENCFF241AOZ_peaks_GRCh38.bed.gz \
| awk '$1=="chr1" && $2<1487904 && $3>1486810' \
> AluY_dup85_400bp_rep1.bed

zcat GSM2423242_ENCFF824IDO_peaks_GRCh38.bed.gz \
| awk '$1=="chr1" && $2<1487904 && $3>1486810' \
> AluY_dup85_400bp_rep2.bed

zcat GSM2423241_ENCFF241AOZ_peaks_GRCh38.bed.gz \
| awk '$1=="chr5" && $2<177489307 && $3>177488206' \
> AluSx1_dup34272_400bp_rep1.bed

zcat GSM2423242_ENCFF824IDO_peaks_GRCh38.bed.gz \
| awk '$1=="chr5" && $2<177489307 && $3>177488206' \
> AluSx1_dup34272_400bp_rep2.bed

zcat GSM2423241_ENCFF241AOZ_peaks_GRCh38.bed.gz \
| awk '$1=="chr5" && $2<177486405 && $3>177485304' \
> AluJr_dup24947_400bp_rep1.bed

zcat GSM2423242_ENCFF824IDO_peaks_GRCh38.bed.gz \
| awk '$1=="chr5" && $2<177486405 && $3>177485304' \
> AluJr_dup24947_400bp_rep2.bed

zcat GSM2423241_ENCFF241AOZ_peaks_GRCh38.bed.gz \
| awk '$1=="chr5" && $2<177466272 && $3>177465174' \
> AluSx_dup34929_400bp_rep1.bed

zcat GSM2423242_ENCFF824IDO_peaks_GRCh38.bed.gz \
| awk '$1=="chr5" && $2<177466272 && $3>177465174' \
> AluSx_dup34929_400bp_rep2.bed

wc -l AluJr_dup74_400bp_rep*.bed \
      AluY_dup85_400bp_rep*.bed \
      AluSx1_dup34272_400bp_rep*.bed \
      AluJr_dup24947_400bp_rep*.bed \
      AluSx_dup34929_400bp_rep*.bed
      
========================================
File: AluJr_dup74_400bp_rep2.bed
chr1	1268613	1268662	HNRNPK_K562_rep02	200	-	1.80861453766413	0.855456925727212	-1	-1
chr1	1268499	1268526	HNRNPK_K562_rep02	200	-	1.80861453766413	0.385857234976744	-1	-1
chr1	1268711	1268727	HNRNPK_K562_rep02	200	-	1.48668644277677	0.183816720322128	-1	-1
========================================
File: AluY_dup85_400bp_rep2.bed
chr1	1487117	1487246	HNRNPK_K562_rep02	200	+	1.57414928402711	1.3426491375667	-1	-1
========================================
File: AluSx1_dup34272_400bp_rep2.bed
chr5	177488999	177489033	HNRNPK_K562_rep02	200	-	2.80861453766413	1.48446100987356	-1	-1
========================================
File: AluSx_dup34929_400bp_rep1.bed
chr5	177466266	177466303	HNRNPK_K562_rep01	200	-	1.04934892080616	0.351195382016362	-1	-1

