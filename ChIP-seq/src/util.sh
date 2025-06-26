#!/bin/bash -O extglob
function run_fastp() {
fastp \
    -i /mnt/e/Alex/X207/RawData/${1}_1.fq.gz \
    -o tmp/${1}_1.fq.gz \
    -I /mnt/e/Alex/X207/RawData/${1}_2.fq.gz \
    -O tmp/${1}_2.fq.gz \
    -w 8 -h Result/06.statistics/fastp_report/${1}.html \
    -j Result/06.statistics/fastp_report/${1}.json
}

function run_bwa() {
    bwa-mem2 mem -t $THREADS $bwa_index tmp/${1}_1.fq.gz tmp/${1}_2.fq.gz | \
    sambamba sort --nthreads 8 tmp/${1}.bam /dev/stdin --out tmp/${1}_sorted.bam
    sambamba markdup --nthreads $THREADS tmp/${1}_sorted.bam | \
    sambamba view --with-header --nthreads $THREADS -f bam \
    --filter "not unmapped and not duplicate and [XA] == null and [XT] == null and mapping_quality >= 20" /dev/stdin Result/02./${1}.bam
}


function run_bowtie2() {
    bowtie2 --threads $THREADS -q --local \
    -U Result/00.raw_data/${1}.fq.gz \
    -x $bowtie2_index | \
    sambamba view --show-progress --sam-input --filter "not unmapped and not duplicate and [XS] == null and mapping_quality >= 30" --nthreads $THREADS --show-progress --format bam /dev/stdin | \
    sambamba sort --show-progress --nthreads $THREADS /dev/stdin --out tmp/${1}.sorted.bam
    sambamba markdup --remove-duplicates --nthreads $THREADS --show-progress tmp/${1}.sorted.bam tmp/${1}.bam 
    bedtools intersect -v -abam tmp/${1}.bam -b $blacklist > Result/02.Alignment/${1}.bam
    sambamba index --show-progress --nthreads $THREADS Result/02.Alignment/${1}.bam
}


function create_bigwig() {
    shopt -s extglob;
    bamCoverage \
    --bam Result/02.Alignment/${1}.bam \
    --outFileName Result/02.Alignment/${1}.bigwig \
    --normalizeUsing RPGC \
    --effectiveGenomeSize $effectiveGenomeSize \
    --extendReads 110 \
    --binSize 20 \
    --blackListFileName $blacklist \
    --numberOfProcessors $THREADS \
    --outFileFormat bigwig
}

# compute the frag length, data quality characteristics based on cross-correlation analysis and/or peak calling
function run_spp() {
    run_spp.R \
        -c=Result/02.Alignment/${1}.bam \
        -odir=Result/04.Visualization/SPP \
        -savp=Result/04.Visualization/SPP/${1}_xcor.pdf \
        -out=Result/04.Visualization/SPP/${1}_xcor_metrics.txt
}

function narrowPeak() {
    local format="${2:-BAM}"
    macs3  callpeak \
    --gsize $effectiveGenomeSize \
    --treatment Result/02.Alignment/${1}.bam \
    --control Result/02.Alignment/${1}.input.bam \
    --bdg --qvalue 0.01 --name ${1} --format $format \
    --outdir Result/03.Peak_Calling/narrowPeak
}
function broadPeak() {
    local format="${2:-BAM}"
    macs3  callpeak \
    --broad --broad-cutoff 0.1 \
    --gsize $effectiveGenomeSize \
    --treatment Result/02.Alignment/${1}.bam \
    --control Result/02.Alignment/${1}.input.bam \
    --bdg --name ${1} --format $format \
    --outdir Result/03.Peak_Calling/broadPeak
}

function run_homer() {
    local format="${2:-narrow}"
    findMotifsGenome.pl Result/03.Peak_Calling/${1}_peaks_${format}Peak_IDR.bed $genome Result/06.Motif/homer/${1}_${format} -size 200 -len 8,10,12
}

function run_meme() {
    local format="${2:-narrow}"

    bedtools getfasta \
    -fi $genome \
    -bed Result/03.Peak_Calling/${1}_peaks_${format}Peak_IDR.bed \
    -fo tmp/${1}_${format}.fna

    meme-chip \
    -oc Result/06.Motif/MEME/${1}_${format} \
    $3 \
    tmp/${1}_${format}.fna
}

