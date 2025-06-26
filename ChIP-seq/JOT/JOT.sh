bwa-mem2 index -p ${outdir}/database/mmu/mmu ${outdir}/database/mmu.v36.gencode.pri.fa.gz
mkdir -p Result/CleanFastq Result/01.fastp
cat RawData/fq.txt | while read i
do
    fastp \
    -i RawData/${i}.fastq.gz \
    -o Result/CleanFastq/${i}.fastq.gz \
    -j Result/01.fastp/${i}.json \
    -h Result/01.fastp/${i}.html
done

cat RawData/fq.txt | while read i;do
bwa-mem2 mem \
    -t 48 \
    ${outdir}/database/mmu/mmu \
    ${outdir}/Result/CleanFastq/${i}.fastq.gz | \
    sambamba view --sam-input --filter "mapping_quality >= 30" --nthreads 48 --show-progress --format bam /dev/stdin --output-filename ${i}.a.bam
    sambamba markdup --remove-duplicates --nthreads 48 --show-progress ${i}.a.bam ${i}.b.bam 
    sambamba sort ${i}.b.bam  --nthreads 48 --show-progress --compression-level 9 --out ${outdir}/Result/02.Align/${i}.bam
    rm ${i}.a.bam ${i}.b.bam
done


# ========= merge all peaks =========
bedops \
--merge ${outdir}/Result/MACS3/hela_1_REST_peaks.narrowPeak ${outdir}/Result/MACS3/hela_2_REST_peaks.narrowPeak \
>${outdir}/Result/MACS3/hela_peaks.bed

# ========= bedgraph to bigwig =========
samtools faidx database/mmu.v36.gencode.pri.fa
for i in WT11 WT12 WT21 WT22 KO1 KO2
do
bedGraphToBigWig Result/MACS3/${i}_control_lambda.bdg  database/mmu.v36.gencode.pri.fa.fai  Result/MACS3/${i}_control.bw
bedGraphToBigWig Result/MACS3/${i}_treat_pileup.bdg  database/mmu.v36.gencode.pri.fa.fai  Result/MACS3/${i}_treat.bw
done