# -*- snakemake -*-
include: "picard.settings.smk"

config_default = {'picard' :{'interval_list_to_bed' : _picard_config_rule_default.copy()}}

update_config(config_default, config)
config = config_default

rule picard_interval_list_to_bed:
    """Picard: convert interval list to bed"""
    params: cmd = config['picard']['cmd'] + INTERVAL_LIST_TO_BED,
            options = config['picard']['interval_list_to_bed']['options'],
            runtime = config['picard']['interval_list_to_bed']['runtime'],
    input: interval_list = "{prefix}.interval_list"
    output: bed = "{prefix}.bed"
    conda: "env.yaml"
    threads: 1
    shell: "{params.cmd} {params.options} I={input.interval_list} O={output.bed}"


localrules: picard_interval_list_to_bed
