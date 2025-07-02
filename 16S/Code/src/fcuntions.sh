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

function print_sample(){
    printf "%s\n" "${samples[@]}"
}