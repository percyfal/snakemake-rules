# -*- snakemake -*-
include: "cutadapt.settings"

rule cutadapt_save_log:
    """Cutadapt: save log output by copying"""
    input: metrics = "{prefix}.cutadapt_metrics"
    output: metrics = "{prefix}.cutadapt_metrics.txt"
    conda: "env.yaml"
    threads: 1
    shell: "cp {input.metrics} {output.metrics}"

localrules: cutadapt_save_log
