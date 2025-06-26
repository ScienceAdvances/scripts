#!/bin/bash
export WORKDIR=/home/zktbk0/Alex/ChIP-seq
cd $WORKDIR
mamba activate chip
source src/util.sh

export DB=/home/zktbk0/Alex/ChIP-seq/DB
export META=meta.tsv
export THREADS=32
export bowtie2_index=/home/zktbk0/Alex/ChIP-seq/DB/mmu
export blacklist=$DB/mm10-blacklist.v2.bed
export effectiveGenomeSize=2654621783
export genome=$DB/mmu.v37.gencode.pri.fna
export gtf=$DB/mmu.v37.gencode.pri.gtf
export -f create_bigwig
export -f run_bowtie2
export -f narrowPeak
export -f run_spp
export -f broadPeak
export -f run_homer
export -f run_meme

mkdir -p logs tmp Result/01.stat/fastp_report Result/02.Alignment Result/03.Peak_Calling/narrowPeak Result/03.Peak_Calling/broadPeak
mkdir -p Result/04.Visualization/SPP  Result/05.Annotation Result/06.Motif/homer Result/06.Motif/MEME

# Get the list of samples from the meta.tsv file
export samples=$(awk 'NR>1 {print $1}' $META)
export non_input_samples=$(grep -v "input" $META | awk 'NR>1 {print $1}')
export groups=$(awk 'NR>1 {if($2!="input")  print $2}' $META  | uniq)
echo ${non_input_samples[@]}

# ========= align to refrence genome =========
printf "%s\n" "${samples[@]}"| parallel -j 1 run_bowtie2

# ========= create bigwig file =========
printf "%s\n" "${samples[@]}"| parallel -j 8 create_bigwig

# ========= SPP =========
printf "%s\n" "${samples[@]}" | parallel  -j 8 run_spp

# call peaks using MACS3
printf "%s\n" "${non_input_samples[@]}" | parallel -j 8 narrowPeak
printf "%s\n" "${non_input_samples[@]}" | parallel -j 8 broadPeak

# ========= consensus peak over replicates =========
# python src/consensus_peak.py
python src/IDR.py

# ========= Annotation =========
Rscript src/05.chipseeker.R

# ========= motif =========
printf "%s\n" "${groups[@]}" | parallel -j 2 run_homer
printf "%s\n" "${groups[@]}" | parallel -j 2 'run_homer {} broad'
printf "%s\n" "${groups[@]}" | parallel -j 2 'run_meme {} narrow "--db $DB/motif/HOCOMOCOv11_full_MOUSE_mono_meme_format.meme --db $DB/motif/JASPAR2024_CORE_vertebrates_non-redundant.meme"'
printf "%s\n" "${groups[@]}" | parallel -j 2 'run_meme {} broad "--db $DB/motif/HOCOMOCOv11_full_MOUSE_mono_meme_format.meme --db $DB/motif/JASPAR2024_CORE_vertebrates_non-redundant.meme"'


# ========= ChIP QC & Diff =========
printf "broad\nnarrow" | parallel -j 2  'Rscript src/07.diffband.R --specie mmu --peak_format {}'


# ========= SPP =========
export -f run_spp
printf "%s\n" "${samples[@]}" | parallel -j 8 run_spp

# >>>>>>>>== ====== deeptools plot ========= >>>>>>>>
# Cumulative enrichment
plotFingerprint \
    --bamfiles $(printf "Result/02.Alignment/%s.bam " ${samples[@]}) \
    --extendReads 110  \
    --binSize 1000 \
    --plotFile Result/04.Visualization/Fingerprint.pdf \
    --outRawCounts Result/04.Visualization/FingerprintRawCounts.tsv \
    --numberOfProcessors $THREADS \
    --labels ${samples[@]}

# Sample clustering
multiBamSummary bins \
    --bamfiles $(printf "Result/02.Alignment/%s.bam " ${samples[@]})  \
    --binSize=5000 \
    --extendReads 110 \
    --minMappingQuality 30 \
    --labels ${samples[@]} \
    --outFileName Result/04.Visualization/CoverageSummary.npz \
    --outRawCounts Result/04.Visualization/CoverageSummary.xls
    # --region 19 \ # limiting the binning of the genome to chromosome 19

# Correlation 
plotCorrelation \
-in Result/04.Visualization/CoverageSummary.npz \
--corMethod pearson --skipZeros \
--plotTitle "Pearson Correlation of Average Scores Per Transcript" \
--whatToPlot scatterplot \
-o Result/04.Visualization/scatter_Pearson_Correlation.pdf   \
--outFileCorMatrix Result/04.Visualization/scatter_Pearson_Correlation.xls

plotCorrelation \
-in Result/04.Visualization/CoverageSummary.npz \
--corMethod pearson --skipZeros \
--plotTitle "Pearson Correlation of Average Scores Per Transcript" \
--whatToPlot heatmap \
-o Result/04.Visualization/heatmap_Pearson_Correlation.pdf   \
--outFileCorMatrix Result/04.Visualization/heatmap_Pearson_Correlation.xls

plotPCA -in Result/04.Visualization/CoverageSummary.npz \
-o Result/04.Visualization/PCA_readCounts.pdf \
-T "PCA of read counts"

# ========= TSS =========
computeMatrix reference-point \
--referencePoint TSS \
--regionsFileName $gtf \
--scoreFileName $(printf "Result/02.Alignment/%s.bigwig " ${samples[@]}) \
--beforeRegionStartLength 1000 \
--afterRegionStartLength 1000 \
--skipZeros \
--outFileName tmp/reference-point-TSS.gz \
--numberOfProcessors $THREADS

plotHeatmap \
--colorMap RdBu \
--matrixFile tmp/reference-point-TSS.gz \
--outFileName Result/04.Visualization/line_heatmap_TSS.pdf \
--whatToShow 'plot, heatmap and colorbar' \
--verbose 

# ========= regions =========
computeMatrix scale-regions \
--regionBodyLength 5000 \
--regionsFileName $gtf \
--scoreFileName $(printf "Result/02.Alignment/%s.bigwig " ${samples[@]}) \
--beforeRegionStartLength 3000 \
--afterRegionStartLength 3000 \
--skipZeros \
--outFileName tmp/scale-regions.gz \
--numberOfProcessors $THREADS

plotHeatmap \
--missingDataColor 1 \
--colorList 'white,#0066CC' \
--heatmapHeight 12 \
--matrixFile tmp/scale-regions.gz \
--outFileName Result/04.Visualization/line_heatmap_regions.pdf \
--whatToShow 'plot, heatmap and colorbar' \
--verbose 

