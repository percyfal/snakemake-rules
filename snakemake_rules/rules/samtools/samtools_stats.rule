# -*- snakemake -*-
include: "samtools.settings"

config_default = {'samtools' :{'stats' : _samtools_config_rule_default.copy()}}

update_config(config_default, config)
config = config_default

rule samtools_stats:
    """BAM statistics"""
    params: options = config['samtools']['stats']['options'],
            cmd = config['samtools']['cmd'],
            runtime = config['samtools']['stats']['runtime']
    input: bam = "{prefix}.bam"
    output: stats = "{prefix}.samtools_stats.txt"
    threads: config['samtools']['stats']['threads']
    conda: "env.yaml"
    shell: "{params.cmd} stats {params.options} {input.bam} > {output.stats}"
