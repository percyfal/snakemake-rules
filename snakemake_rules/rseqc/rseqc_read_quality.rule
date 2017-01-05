# -*- snakemake -*-
include: "rseqc.settings"

config_default = {'rseqc' :{'read_quality' : _rseqc_config_rule_default.copy()}}

update_config(config_default, config)
config = config_default


rule rseqc_read_quality:
    """Run RSeQC read_quality.py"""
    params: cmd = READ_QUALITY,
            options = config['rseqc']['read_quality']['options'],
            runtime = config['rseqc']['read_quality']['runtime'],
    input: "{prefix}.bam"
    output: r = "{prefix}_rseqc/read_quality.qual.r"
    threads: config['rseqc']['read_quality']['threads']
    conda: "env.yaml"
    shell: "{params.cmd} {params.options} -i {input} -o $(dirname {output.r})/read_quality"

