#/usr/bin/env bash
export WORKDIR=/home/alex/MetaFlow
cd $WORKDIR
source $WORKDIR/scripts/function.sh
export -f $(declare -F | awk '{print $NF}')
export MetaFlowDB=/mnt/e/Alex/X207/MetaFlowDB
META=/mnt/e/Alex/X207/meta.tsv
export THREADS=45
export bwa_index=/mnt/e/Alex/DATAHUB/index/hsa_bwa/hsa
export bowtie2_index=/mnt/e/Alex/DATAHUB/index/hsa_bowtie2/hsa
megan=/home/data/wd/Vik/APP/megan/tools
mamba activate MetaFlow



mkdir -p tmp Result/01.clean_data/ Result/02.Alignment tmp Result/03.Peak_Calling/NarrowPeak Result/03.Peak_Calling/BroadPeak Result/04.Visualization Result/05.Annotation Result/10.stat/fastp_report

# Get the list of samples from the meta.tsv file
mapfile -t samples < <(awk 'NR>1 {print $1}' $META)
export samples

# Run fastp on each sample
print_sample | parallel -j 16 run_fastp

# bwa-mem2: map reads to genome, sambamba markdup: remove PCR duplicates, sambamba vie: sort and filter
print_sample | parallel -j 4 run_bwa

# remove black list
print_sample | parallel -j 16 "bedtools intersect -v -abam tmp/{}.bam -b src/mm10-blacklist.v2.bed > Result/02.Alignment/{}.bam"

# call peaks using MACS3
macs3 callpeak -g hs --treatment Result/02.Alignment/{}.bam --control Result/02.Alignment/{}.bam --pvalue 0.05 --name {} --format BAMPE --outdir Result/03.Peak_Calling/NarrowPeak
macs3 callpeak -g hs --treatment Result/02.Alignment/{}.bam --control Result/02.Alignment/{}.bam --broad --broad-cutoff 0.1 --name {} --format BAMPE --outdir Result/03.Peak_Calling/BroadPeak


# 
print_sample | parallel -j 4 bamCoverage -b Result/02.Alignment/{}.bam -o cResult/02.Alignment/{}.bw



# homer Annotaion


# stats
multiBamSummary bins \
    --bamfiles $(printf "Result/02.Alignment/*bam " "${samples[@]}")  \
    --minMappingQuality 30 \
    --labels ${samples[@]} \
    -out Result/10.stat/CoverageSummary.npz \
    --outRawCounts Result/10.stat/CoverageSummary.xls
    # --region 19 \ # limiting the binning of the genome to chromosome 19

# Calculate FRiP score https://deeptools.readthedocs.io/en/develop/content/example_api_tutorial.html#computing-the-frip-score
# 一个高质量的文库其FRiP score值应该大于0.03，最低也要大于0.02
import deeptools
import pysam
bed_file = "Result/03.Peak_Calling/NarrowPeak/{}_peaks.narrowPeak"
bam_file = "Result/02.Alignment/{}.bam"
cr = deeptools.countReadsPerBin.CountReadsPerBin(bam_file, bedFile=bed_file, numberOfProcessors=10)
reads_at_peaks = cr.run()
total = reads_at_peaks.sum(axis=0)
bam = pysam.AlignmentFile(bam_file)
frip = float(total) / bam.mapped
pd.DataFrame(dict(Sample=samples,FRiP=frips)).to_excel("Result/10.stat/FRiP_score.xlsx",index=False)

# Correlation 
plotCorrelation \
-in Result/10.stat/CoverageSummary.npz \
--corMethod pearson --skipZeros \
--plotTitle "Pearson Correlation of Average Scores Per Transcript" \
--whatToPlot scatterplot \
-o Result/10.stat/scatter_Pearson_Correlation.pdf   \
--outFileCorMatrix Result/10.stat/scatter_Pearson_Correlation.xls

plotCorrelation \
-in Result/10.stat/CoverageSummary.npz \
--corMethod pearson --skipZeros \
--plotTitle "Pearson Correlation of Average Scores Per Transcript" \
--whatToPlot heatmap \
-o Result/10.stat/heatmap_Pearson_Correlation.pdf   \
--outFileCorMatrix Result/10.stat/heatmap_Pearson_Correlation.xls

plotPCA -in Result/10.stat/CoverageSummary.npz \
-o Result/10.stat/PCA_readCounts.pdf \
-T "PCA of read counts"


idr \
--samples $(printf "Result/03.Peak_Calling/NarrowPeak/%s_peaks.narrowPeak","${samples[@]}") \
--input-file-type narrowPeak \
--rank p.value \
--output-file IDR_narrowPeak \
--plot \
--log-output-file IDR_narrowPeak.log

idr \
--samples $(printf "Result/03.Peak_Calling/NarrowPeak/%s_peaks.broadPeak","${samples[@]}") \
--input-file-type broadPeak \
--rank p.value \
--output-file IDR_broadPeak \
--plot \
--log-output-file IDR_broadPeak.log
