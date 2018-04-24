# -*- snakemake -*-
include: "samtools.settings.smk"

config_default = {'samtools' :{'faidx' : _samtools_config_rule_default.copy()}}

update_config(config_default, config)
config = config_default


rule samtools_faidx:
    """Run samtools faidx"""
    params: cmd = config['samtools']['cmd'],
            options = config['samtools']['faidx']['options'],
            runtime = config['samtools']['faidx']['runtime']
    input: "{prefix}.fa"
    output: "{prefix}.fa.fai"
    threads: config['samtools']['faidx']['threads']
    conda: "env.yaml"
    shell: "{params.cmd} {params.options} faidx {input}"
