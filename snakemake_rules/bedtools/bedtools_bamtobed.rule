# -*- snakemake -*-
include: "bedtools.settings"

config_default = {'bedtools' :{'bamtobed' : _bedtools_config_rule_default.copy()}}

update_config(config_default, config)
config = config_default


rule bedtools_bamtobed:
    """Convert a bam file to bed format"""
    params: cmd = " ".join([BEDTOOLS, BEDTOOLS_BAMTOBED]),
            options = config['bedtools']['bamtobed']['options'],
            runtime = config['bedtools']['bamtobed']['runtime'],
    input: "{prefix}.bam"
    output: "{prefix}.bed"
    threads: config['bedtools']['bamtobed']['threads']
    conda: "env.yaml"
    shell: "{params.cmd} {params.options} -i {input} > {output}"
