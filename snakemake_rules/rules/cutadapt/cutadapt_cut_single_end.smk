# -*- snakemake -*-
include: "cutadapt.settings.smk"

config_default = {'cutadapt' :{'single_end' : _cutadapt_config_rule_default.copy()}}

update_config(config_default, config)
config = config_default

rule cutadapt_cut_single_end:
    """Cutadapt: cut single end sequences. 

    Trim both fiveprime and threeprime adatper.
    """
    params: cmd = config['cutadapt']['cmd'],
            options = config['cutadapt']['single_end']['options'],
            threeprime = config['cutadapt']['threeprime'],
            fiveprime = config['cutadapt']['fiveprime'],
            runtime = config['cutadapt']['runtime']
    wildcard_constraints: fastq = config['cutadapt']['fastq_suffix_re'],
                          read = config['cutadapt']['read_label_re']
    input: read = "{prefix}{read}{fastq}"
    output: read = "{prefix}.trimmed{read}{fastq}",
            log = "{prefix}.trimmed{read}{fastq}.cutadapt_metrics"
    threads: config['cutadapt']['threads']
    conda: "env.yaml"
    log: "{prefix}.trimmed{read}{fastq}.cutadapt_metrics"
    shell: "{params.cmd} {params.options} {input.read} -a {params.threeprime} -o {output.read} > {log}"
