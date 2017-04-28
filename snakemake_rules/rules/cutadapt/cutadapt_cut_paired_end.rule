# -*- snakemake -*-
include: "cutadapt.settings"

config_default = {'cutadapt' :{'paired_end' : _cutadapt_config_rule_default.copy()}}

update_config(config_default, config)
config = config_default


rule cutadapt_cut_paired_end:
    """Cutadapt: cut paired end sequences
    
    Cut both five- and threeprime adapter.
    """
    params: cmd = config['cutadapt']['cmd'],
            options = config['cutadapt']['paired_end']['options'],
            threeprime = config['cutadapt']['threeprime'],
            fiveprime = config['cutadapt']['fiveprime'],
            runtime = config['cutadapt']['runtime'],
    wildcard_constraints: fastq = config['cutadapt']['fastq_suffix_re']
    input: read1 = "{prefix}" + config['cutadapt']['read1_label'] + "{fastq}",
           read2 = "{prefix}" + config['cutadapt']['read2_label'] + "{fastq}"
    output:
        read1 = "{prefix}.trimmed" + config['cutadapt']['read1_label'] + "{fastq}",
        read2 = "{prefix}.trimmed" + config['cutadapt']['read2_label'] + "{fastq}",
        log = "{prefix}{fastq}.cutadapt_metrics"
    conda: "env.yaml"
    threads: config['cutadapt']['threads']
    log: "{prefix}{fastq}.cutadapt_metrics"
    shell: "{params.cmd} {params.options} {input.read1} {input.read2} -a {params.threeprime} -A {params.fiveprime} -o {output.read1} -p {output.read2} > {log}"
