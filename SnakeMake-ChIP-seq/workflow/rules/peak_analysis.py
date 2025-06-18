rule macs_callpeak_broad:
    input:
        treatment="results/filtered/{sample}.sorted.bam",
        control="results/filtered/{control}.sorted.bam",
        gsizepath="resources/ref/gsize.txt"
    output:
        # all output-files must share the same basename and only differ by it's extension
        # Usable extensions (and which tools they implicitly call) are listed here:
        #         https://snakemake-wrappers.readthedocs.io/en/stable/wrappers/macs2/callpeak.html.
        multiext("results/macs_callpeak/{sample}-{control}.broad",
                 "_peaks.xls",
                 # these output extensions internally set the --bdg or -B option:
                 "_treat_pileup.bdg",
                 "_control_lambda.bdg",
                 # these output extensions internally set the --broad option:
                 "_peaks.broadPeak",
                 "_peaks.gappedPeak"
                 )
    log:
        "logs/macs/callpeak.{sample}-{control}.broad.log"
    params:
        name = "{sample}",
        format = "BAMPE" if config["is_pe"] else "BAM",
        pvalue = config["params"]['callpeak']['pvalue'],
        qvalue = config["params"]['callpeak']['qvalue'],
        keep_dup = config["params"]['callpeak']['keep-dup'],
        nolambda = config["params"]['callpeak']['nolambda'],
        broad = config["params"]['callpeak']['broad'],
        broad_cutoff = config["params"]['callpeak']['broad-cutoff'],
    conda: "../envs/call_peak.yml"
    scripts: "../scripts/call_peak.py"

rule macs_callpeak_narrow:
    input:
        treatment="results/filtered/{sample}.sorted.bam",
        control="results/filtered/{control}.sorted.bam",
        gsizepath="resources/ref/gsize.txt"
    output:
        # all output-files must share the same basename and only differ by it's extension
        # Usable extensions (and which tools they implicitly call) are listed here:
        #         https://snakemake-wrappers.readthedocs.io/en/stable/wrappers/macs2/callpeak.html.
        multiext("results/macs_callpeak/{sample}-{control}.narrow",
                 "_peaks.xls",
                 # these output extensions internally set the --bdg or -B option:
                 "_treat_pileup.bdg",
                 "_control_lambda.bdg",
                 # these output extensions internally set the --broad option:
                 "_peaks.narrowPeak",
                 "_summits.bed"
                 )
    log:
        "logs/macs/callpeak.{sample}-{control}.narrow.log"
    params:
        lambda w, input: "-f {bam_format} {gsize} -B --SPMR --keep-dup all {pvalue} {qvalue}".format(
            gsize="{}".format(open(input.gsizepath).read().strip()),
            pvalue="-p {}".format(config["params"]["callpeak"]["p-value"]) if config["params"]["callpeak"]["p-value"] else "",
            qvalue="-q {}".format(config["params"]["callpeak"]["q-value"]) if config["params"]["callpeak"]["q-value"] else "",
            bam_format="BAM" if config["single_end"] else "BAMPE")
    conda: "../envs/call_peak.yml"
    scripts: "../scripts/call_peak.py"

rule peaks_count:
    input:
        peaks="results/macs_callpeak/{sample}-{control}.{peak}_peaks.{peak}Peak"
    output:
        "results/macs_callpeak/peaks_count/{sample}-{control}.{peak}.peaks_count.tsv"
    log:
        "logs/macs_callpeak/peaks_count/{sample}-{control}.{peak}.peaks_count.log"
    conda:
        "../envs/gawk.yaml"
    shell:
        "cat {input.peaks} | "
        " wc -l | "
        " gawk -v OFS='\t' '{{print \"{wildcards.sample}-{wildcards.control}_{wildcards.peak}_peaks\", $1}}' "
        " > {output} 2> {log}"

rule sm_report_peaks_count_plot:
    input:
        get_peaks_count_plot_input()
    output:
        report("results/macs_callpeak/plots/plot_{peak}_peaks_count.pdf", caption="../report/plot_peaks_count_macs.rst", category="CallPeaks")
    log:
        "logs/macs_callpeak/plot_{peak}_peaks_count.log"
    conda:
        "../envs/r_plots.yaml"
    script:
        "../scripts/plot_peaks_count_macs.R"

rule bedtools_intersect:
    input:
        left="results/filtered/{sample}.sorted.bam",
        right="results/macs_callpeak/{sample}-{control}.{peak}_peaks.{peak}Peak"
    output:
        pipe("results/bedtools_intersect/{sample}-{control}.{peak}.intersected.bed")
    params:
        extra="-bed -c -f 0.20"
    log:
        "logs/bedtools/intersect/{sample}-{control}.{peak}.intersected.log"
    wrapper:
        config["warpper_mirror"] + "bio/bedtools/intersect"

rule frip_score:
    input:
        intersect="results/bedtools_intersect/{sample}-{control}.{peak}.intersected.bed",
        flagstats=expand("results/{step}/{{sample}}.sorted.{step}.flagstat", step= "bamtools_filtered" if config["single_end"]
        else "orph_rm_pe")
    output:
        "results/bedtools_intersect/{sample}-{control}.{peak}.peaks_frip.tsv"
    log:
        "logs/bedtools/intersect/{sample}-{control}.{peak}.peaks_frip.log"
    conda:
        "../envs/gawk.yaml"
    shell:
        "grep 'mapped (' {input.flagstats} | "
        " gawk -v a=$(gawk -F '\t' '{{sum += $NF}} END {{print sum}}' < {input.intersect}) "
        " -v OFS='\t' "
        " '{{print \"{wildcards.sample}-{wildcards.control}_{wildcards.peak}_peaks\", a/$1}}' "
        " > {output} 2> {log}"

rule sm_rep_frip_score:
    input:
        get_frip_score_input()
    output:
        report("results/macs_callpeak/plots/plot_{peak}_peaks_frip_score.pdf", caption="../report/plot_frip_score_macs_bedtools.rst", category="CallPeaks")
    log:
        "logs/bedtools/intersect/plot_{peak}_peaks_frip_score.log"
    conda:
        "../envs/r_plots.yaml"
    script:
        "../scripts/plot_frip_score.R"

rule create_igv_peaks:
    input:
        "results/macs_callpeak/{sample}-{control}.{peak}_peaks.{peak}Peak"
    output:
        "results/IGV/macs_callpeak-{peak}/merged_library.{sample}-{control}.{peak}_peaks.igv.txt"
    log:
        "logs/igv/create_igv_peaks/merged_library.{sample}-{control}.{peak}_peaks.log"
    shell:
        " find {input} -type f -name '*_peaks.{wildcards.peak}Peak' -exec echo -e 'results/IGV/macs_callpeak/{wildcards.peak}/\"{{}}\"\t0,0,178' \; > {output} 2> {log}"

rule homer_annotatepeaks:
    input:
        peaks="results/macs_callpeak/{sample}-{control}.{peak}_peaks.{peak}Peak",
        genome="resources/ref/genome.fasta",
        gtf="resources/ref/annotation.gtf"
    output:
        annotations="results/homer/annotate_peaks/{sample}-{control}.{peak}_peaks.annotatePeaks.txt"
    threads: 2
    params:
        mode="",
        extra="-gid"
    log:
        "logs/homer/annotate_peaks/{sample}-{control}.{peak}.log"
    wrapper:
        config["warpper_mirror"] + "bio/homer/annotatePeaks"

rule plot_macs_qc:
    input:
        get_macs_peaks()
    output:  #ToDo: add description to report caption
        summmary="results/macs_callpeak/plots/plot_{peak}_peaks_macs_summary.txt",
        plot=report("results/macs_callpeak/plots/plot_{peak}_peaks_macs.pdf", caption="../report/plot_macs_qc.rst", category="CallPeaks")
    params:
        input = lambda wc, input: ','.join(input),
        sample_control_combinations = ','.join(get_sample_control_peak_combinations_list())
    log:
        "logs/macs_callpeak/plot_{peak}_peaks_macs.log"
    conda:
        "../envs/plot_macs_annot.yaml"
    shell:
        "Rscript ../workflow/scripts/plot_macs_qc.R -i {params.input} -s {params.sample_control_combinations}  -o {output.plot} -p {output.summmary} 2> {log}"

rule plot_homer_annotatepeaks:
    input:
        get_plot_homer_annotatepeaks_input()
    output:  #ToDo: add description to report caption
        summmary="results/homer/plots/plot_{peak}_annotatepeaks_summary.txt",
        plot=report("results/homer/plots/plot_{peak}_annotatepeaks.pdf", caption="../report/plot_annotatepeaks_homer.rst", category="CallPeaks")
    params:
        input = lambda wc, input: ','.join(input),
        sample_control_combinations = ','.join(get_sample_control_peak_combinations_list())
    log:
        "logs/homer/plot_{peak}_annotatepeaks.log"
    conda:
        "../envs/plot_macs_annot.yaml"
    shell:
        "Rscript ../workflow/scripts/plot_homer_annotatepeaks.R -i {params.input} -s {params.sample_control_combinations}  -o {output.plot} -p {output.summmary} 2> {log}"

rule plot_sum_annotatepeaks:
    input:
        "results/homer/plots/plot_{peak}_annotatepeaks_summary.txt"
    output:
        report("results/homer/plots/plot_{peak}_annotatepeaks_summary.pdf", caption="../report/plot_annotatepeaks_summary_homer.rst", category="CallPeaks")
    log:
        "logs/homer/plot_{peak}_annotatepeaks_summary.log"
    conda:
        "../envs/r_plots.yaml"
    script:
        "../scripts/plot_annotatepeaks_summary_homer.R"
