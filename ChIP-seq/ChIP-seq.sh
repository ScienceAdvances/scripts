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



mkdir -p tmp Result/01.clean_data/ Result/02.Alignment tmp Result/03.Peak_Calling/NarrowPeak Result/03.Peak_Calling/BroadPeak Result/04.Visualization Result/05.Annotation Result/10.statistics/fastp_report

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

