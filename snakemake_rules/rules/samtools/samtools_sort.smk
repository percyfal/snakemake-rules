# -*- snakemake -*-
include: "samtools.settings"

config_default = {'samtools' :{'sort' : _samtools_config_rule_default.copy()}}

update_config(config_default, config)
config = config_default


rule samtools_sort:
    """Run samtools sort"""
    params: options = config['samtools']['sort']['options'],
            cmd = config['samtools']['cmd'],
            runtime = config['samtools']['sort']['runtime']
    input: "{prefix}.bam"
    output: "{prefix}.sort.bam"
    threads: config['samtools']['sort']['threads']
    conda: "env.yaml"
    shell: "{params.cmd} sort {params.options} -T /tmp/$(basename {wildcards.prefix}).sorted -o {output} {input}"
