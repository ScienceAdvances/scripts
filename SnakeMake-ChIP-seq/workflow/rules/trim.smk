if config["fastqs"].get("pe"):
    rule fastp_pe:
        input:
            sample=get_fastq
        output:
            trimmed=[temp("results/trimmed/{s}_{u}_1.fastq.gz"), temp("results/trimmed/{s}_{u}_2.fastq.gz")],
            html=temp("report/{s}_{u}.fastp.html"),
            json=temp("report/{s}_{u}.fastp.json"),
        log:
            "logs/fastp/{s}_{u}.log"
        threads: 16
        wrapper:
            config["warpper_mirror"]+"bio/fastp"
else:
    rule fastp_se:
        input:
            sample=get_fastq
        output:
            trimmed=temp("results/trimmed/{s}_{u}.fastq.gz"),
            html=temp("report/{s}_{u}.fastp.html"),
            json=temp("report/{s}_{u}.fastp.json"),
        log:
            "logs/fastp/{s}_{u}.log"
        threads: 16
        wrapper:
            config["warpper_mirror"]+"bio/fastp"