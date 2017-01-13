# -*- snakemake -*-
include: "samtools.settings"

config_default = {'samtools' :{'depth' : _samtools_config_rule_default.copy()}}

update_config(config_default, config)
config = config_default


rule samtools_depth:
    """Run samtools depth"""
    params: cmd = config['samtools']['cmd'],
            options = config['samtools']['depth']['options'],
            runtime = config['samtools']['depth']['runtime']
    input: fofn = "{prefix}.fofn"
    output: coverage = "{prefix}.coverage.txt"
    threads: config['samtools']['depth']['threads']
    conda: "env.yaml"
    shell: "{params.cmd} {params.options} depth -f {input.fofn} > {output.coverage}"
