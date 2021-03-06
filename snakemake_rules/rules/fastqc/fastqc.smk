# -*- snakemake -*-
include: "fastqc.settings.smk"


rule fastqc:
    """fastqc: run fastqc on a fastq file"""
    params: cmd = config['fastqc']['cmd'],
            options = config['fastqc']['options'],
            runtime = config['fastqc']['runtime']
    threads: config['fastqc']['threads']
    input: fastq = "{prefix}" + config["ngs.settings"]["fastq_suffix"]
    output: html = "{prefix}_fastqc/fastqc_report.html",
            txt="{prefix}_fastqc/fastqc_data.txt",
            zipfile="{prefix}_fastqc.zip"
    conda: "env.yaml"
    shell:
        "mkdir -p $(dirname {output.html}) && " + \
            "{params.cmd} -t {threads} {params.options} --extract {input} -o $(dirname $(dirname {output.html}))"

