rule test_bowtie2:
    input:
        sample=["reads/{sample}.1.fastq", "reads/{sample}.2.fastq"],
        ref="genome.fasta", #Required for CRAM output
        idx=multiext(
            "index/genome",
            ".1.bt2",
            ".2.bt2",
            ".3.bt2",
            ".4.bt2",
            ".rev.1.bt2",
            ".rev.2.bt2",
        ),

    output:
        "results/mapped/{sample}.cram",
    i   dx="mapped_idx/{sample}.cram.bai",
        metrics="mapped_idx/{sample}.metrics.txt",
        unaligned="mapped_idx/{sample}.unaligned.sam",
        unpaired="mapped_idx/{sample}.unpaired.sam",
        # unconcordant="",
        # concordant="",

    log:
        "logs/bowtie2/{sample}.log",
    params:
        extra="",  # optional parameters
    threads: 8  # Use at least two threads
    wrapper:
        config["warpper_mirror"] + "bio/bowtie2/align"


rule bwa_mem:
    input:
        reads = get_map_reads_input,
        idx = rules.bwa_index.output
    output:
        temp("results/mapped/{sample}-{unit}.bam")
    log:
        "logs/bwa/bwa_mem/{sample}-{unit}.log"
    params:
        index= lambda w, input: os.path.splitext(input.idx[0])[0],
        extra= get_read_group,
        sort="samtools",
        sort_order="coordinate",
    threads: 8
    wrapper:
        config["warpper_mirror"] + "bio/bwa/mem"

rule merge_bams:
    input:
        lambda w: expand("results/mapped/{sample}-{unit}.bam",
            sample = w.sample,
            unit = units.loc[units['sample'] == w.sample].unit.to_list()
        )
    output:
        temp("results/merged/{sample}.bam")
    log:
        "logs/picard/mergebamfiles/{sample}.log"
    params:
        "VALIDATION_STRINGENCY=LENIENT SORT_ORDER=coordinate"
    wrapper:
        config["warpper_mirror"] + "bio/picard/mergesamfiles"

rule mark_merged_duplicates:
    input:
        "results/merged/{sample}.bam"
    output:
        bam=temp("results/picard_dedup/{sample}.bam"),
        metrics="results/picard_dedup/{sample}.metrics.txt"
    log:
        "logs/picard/picard_dedup/{sample}.log"
    params:
        "REMOVE_DUPLICATES=false ASSUME_SORTED=true VALIDATION_STRINGENCY=LENIENT"
    wrapper:
        config["warpper_mirror"] + "bio/picard/markduplicates"
