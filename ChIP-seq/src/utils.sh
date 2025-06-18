#/usr/bin/env bash
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
