# -*- snakemake -*-
include: "cutadapt.settings"

rule cutadapt_cut_fiveprime:
    """Cutadapt: cut fiveprime adapter"""
    params: cmd = config['cutadapt']['cmd'],
            options = config['cutadapt']['options'],
            fiveprime = config['cutadapt']['fiveprime'],
            runtime = config['cutadapt']['runtime'],
    wildcard_constraints: read = config['cutadapt']['read_label_re'],
                          fastq = config['cutadapt']['fastq_suffix_re']
    input: read = "{prefix}{read}{fastq}"
    output: read = "{prefix}.trimmed{read}{fastq}",
            log = "{prefix}.trimmed{read}{fastq}.cutadapt_metrics"
    conda: "env.yaml"
    threads: config['cutadapt']['threads']
    log: "{prefix}.trimmed{read}{fastq}.cutadapt_metrics"
    shell: "{params.cmd} {params.options} -a {params.fiveprime} {input.read} -o {output.read} > {log}"

