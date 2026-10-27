pyGenomeTracks \
--tracks tracks.ini \
--region chr12:47775000-47840000 \
--dpi 300 \
--outFileName HDAC7_eCLIP_tracks.pdf


pyGenomeTracks     --tracks tracks.ini     --region chr12:47775000-47840000     --dpi 300     --outFileName HDAC7_eCLIP_bigWig_all_peaks_2.pdf


 1000  pyGenomeTracks     --tracks tracks_merge.ini     --region chr12:47775000-47840000     --dpi 300     --outFileName HDAC7_eCLIP_bigWig_merged.pdf
 1001  cd projects/ASE_NB_project/eCLIP/scripts/
 1002  pyGenomeTracks     --tracks tracks_merge.ini     --region chr12:47775000-47840000     --dpi 300     --outFileName HDAC7_eCLIP_bigWig_merged.pdf
 1003  pyGenomeTracks     --tracks tracks_merge.ini     --region chr12:47775000-47840000     --dpi 300     --outFileName HDAC7_eCLIP_bigWig_merged_filtered.pdf
 1004  cat /home/2024rbb1003/projects/ASE_NB_project/eCLIP/results/HDAC7_peaks_p0.05_rep1.bed /home/2024rbb1003/projects/ASE_NB_project/eCLIP/results/HDAC7_peaks_p0.05_rep2.bed > /home/2024rbb1003/projects/ASE_NB_project/eCLIP/results/HDAC7_combined.bed
 1005  sort -k1,1 -k2,2n /home/2024rbb1003/projects/ASE_NB_project/eCLIP/results/HDAC7_combined.bed > /home/2024rbb1003/projects/ASE_NB_project/eCLIP/results/HDAC7_combined.sorted.bed
 1006  bedtools merge -i /home/2024rbb1003/projects/ASE_NB_project/eCLIP/results/HDAC7_combined.sorted.bed > /home/2024rbb1003/projects/ASE_NB_project/eCLIP/results/HDAC7_consensus_peaks.bed
 1007  bedtools merge -i HDAC7_combined.sorted.bed -c 4,5,6 -o collapse,max,distinct > HDAC7_consensus_peaks_annotated.bed
 1008  bedtools merge -i /home/2024rbb1003/projects/ASE_NB_project/eCLIP/results/HDAC7_combined.sorted.bed -c 4,5,6 -o collapse,max,distinct > HDAC7_consensus_peaks_annotated.bed
 1009  cd ..
 1010  cd results/
 1011  bedtools intersect -a HDAC7_peaks_p0.05_rep1.bed -b HDAC7_peaks_p0.05_rep2.bed > HDAC7_shared_peaks.bed

conda activate pygenome_env
pyGenomeTracks     --tracks tracks_consensus.ini     --region chr12:47775000-47840000     --dpi 300     --outFileName HDAC7_eCLIP_bigWig_consensus.pdf
pyGenomeTracks     --tracks tracks_consensus.ini     --region chr12:47775000-47840000     --dpi 300     --outFileName HDAC7_eCLIP_bigWig_shared_2.pdf
pyGenomeTracks     --tracks tracks_consensus.ini     --region chr12:47790000-47800000     --dpi 300     --outFileName HDAC7_eCLIP_bigWig_shared_3.pdf
pyGenomeTracks     --tracks tracks_consensus.ini     --region chr12:47789812-47791308     --dpi 300     --outFileName HDAC7_eCLIP_bigWig_shared_zoom.pdf
pyGenomeTracks     --tracks tracks_consensus.ini     --region chr12:47780000-47800000     --dpi 300     --outFileName HDAC7_eCLIP_bigWig_shared_4.pdf
pyGenomeTracks     --tracks tracks_consensus.ini     --region chr12:47785000-47800000     --dpi 300     --outFileName HDAC7_eCLIP_bigWig_shared_5.pdf
pyGenomeTracks     --tracks tracks_consensus.ini     --region chr12:47789000-47797000     --dpi 300     --outFileName HDAC7_eCLIP_bigWig_shared_5.pdf
pyGenomeTracks     --tracks tracks_consensus.ini     --region chr12:47790000-47797000     --dpi 300     --outFileName HDAC7_eCLIP_bigWig_shared_6.pdf


####Alu
pyGenomeTracks     --tracks tracks_AluY_dup85_rep2.ini     --region chr1:1486810-1487904     --dpi 300     --outFileName AluY_dup85_rep2_eCLIP_bigWig.pdf

pyGenomeTracks     --tracks tracks_AluJr_dup74_400bp_rep2.ini     --region chr1:1268296-1269365    --dpi 300     --outFileName AluJr_dup74_400bp_rep2_eCLIP_bigWig.pdf

pyGenomeTracks     --tracks tracks_AluSx1_dup34272_400bp_rep2.ini     --region chr5:177488206-177489307    --dpi 300     --outFileName AluSx1_dup34272_400bp_rep2_eCLIP_bigWig.pdf

pyGenomeTracks     --tracks tracks_AluSx_dup34929_400bp_rep1.ini     --region chr5:177465174-177466272    --dpi 300     --outFileName AluSx_dup34929_400bp_rep1_eCLIP_bigWig.pdf
