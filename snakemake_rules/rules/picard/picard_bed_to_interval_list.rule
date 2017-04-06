# -*- snakemake -*-
include: "picard.settings"
include: "picard_create_sequence_dictionary.rule"

config_default = {'picard' :{'bed_to_interval_list' : _picard_config_rule_default.copy()}}

update_config(config_default, config)
config = config_default

rule picard_bed_to_interval_list:
    """Picard: convert bed file to interval list"""
    params: cmd = config['picard']['cmd'] + BED_TO_INTERVAL_LIST,
            options = config['picard']['bed_to_interval_list']['options'],
            runtime = config['picard']['bed_to_interval_list']['runtime']
    input: bed = "{prefix}.bed", dict = os.path.splitext(config['picard']['ref'])[0] + ".dict"
    output: list = "{prefix}.interval_list"
    conda: "env.yaml"
    threads: config['picard']['bed_to_interval_list']['threads']
    shell: "{params.cmd} {params.options} I={input.bed} O={output.list} SD={input.dict}"

