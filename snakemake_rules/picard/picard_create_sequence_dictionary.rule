# -*- snakemake -*-
include: "picard.settings"

config_default = {'picard' :{'create_sequence_dictionary' : _picard_config_rule_default.copy()}}

update_config(config_default, config)
config = config_default


rule picard_create_sequence_dictionary:
    """Picard: create sequence dictionary"""
    params: cmd = config['picard']['cmd'] + CREATE_SEQUENCE_DICTIONARY,
            options = config['picard']['create_sequence_dictionary']['options'],
            runtime = config['picard']['create_sequence_dictionary']['runtime']
    input: fa="{prefix}.fa"
    output: dict="{prefix}.dict"
    conda: "env.yaml"
    threads: config['picard']['create_sequence_dictionary']['threads']
    shell: "{params.cmd} {params.options} R={input.fa} O={output.dict}"

