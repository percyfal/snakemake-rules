# -*- snakemake -*-
include: "picard.settings.smk"
include: "picard_build_bam_index.smk"

config_default = {'picard' :{'collect_insert_size_metrics' : _picard_config_rule_default.copy()}}
config_default['picard']['collect_insert_size_metrics'].update({'targets' : []})

update_config(config_default, config)
config = config_default


rule picard_collect_insert_size_metrics:
    """Picard: collect insertion size metrics"""
    params: cmd = config['picard']['cmd'] + COLLECT_INSERT_SIZE_METRICS,
            options = config['picard']['collect_insert_size_metrics']['options'],
            suffix = ".insert_metrics",
            runtime = config['picard']['collect_insert_size_metrics']['runtime']
    input: bam = "{prefix}.bam", bai = "{prefix}.bai", ref = config['picard']['ref']
    output: metrics = "{prefix}.insert_metrics"
    conda: "env.yaml"
    threads: config['picard']['collect_insert_size_metrics']['threads']
    shell: "{params.cmd} {params.options} H={wildcards.prefix}.hist I={input.bam} O={output.metrics} R={input.ref}"

