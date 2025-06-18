rule samtools_index:
    input:
        "results/{step}/{samples_units}.bam"
    output:
        "results/{step}/{samples_units}.bam.bai"
    params:
        "" # optional params string
    log:
        "logs/samtools-index/{step}/{samples_units}.log"
    wrapper:
        config["warpper_mirror"] + "bio/samtools/index"
