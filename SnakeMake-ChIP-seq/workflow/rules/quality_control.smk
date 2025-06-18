rule plot_fingerprint:
    input:
        bam_files=["results/filtered/{sample}.sorted.bam", "results/filtered/{control}.sorted.bam"],
        bam_idx=["results/filtered/{sample}.sorted.bam.bai", "results/filtered/{control}.sorted.bam.bai"],
        jsd_sample="results/filtered/{control}.sorted.bam",
        stats=expand("results/{step}/{{sample}}.sorted.{step}.stats.txt", step="orph_rm_pe" if config["is_pe"] else "bamtools_filtered")
    output:
        fingerprint=report("results/deeptools/{sample}-{control}.plot_fingerprint.pdf", caption="../report/plot_fingerprint_deeptools.rst", category="QC"),
        counts="results/deeptools/{sample}-{control}.fingerprint_counts.txt",
        qc_metrics="results/deeptools/{sample}-{control}.fingerprint_qcmetrics.txt"
    log:
        "logs/deeptools/plot_fingerprint.{sample}-{control}.log"
    params:
        "--labels {sample} {control}",
        "--skipZeros ",
        "--numberOfSamples 500000.0 ",
        lambda w, input:
            "{se_option}{fragment_size}".format(
                se_option="--extendReads " if not config["is_pe"] else "",
                # Estimated fragment size used to extend single-end reads
                fragment_size=
                "$(grep ^SN {stats} | "
                "cut -f 2- | "
                "grep -m1 'average length:' | "
                "awk '{{print $NF}}') ".format(
                    stats=input.stats)
                if not config["is_pe"] else ""
            )
    threads: 8
    wrapper:
        config["warpper_mirror"] + "bio/deeptools/plotfingerprint"