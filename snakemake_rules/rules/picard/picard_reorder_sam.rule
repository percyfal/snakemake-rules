# -*- snakemake -*-
include: "picard.settings"

config_default = {'picard' :{'reorder_sam' : _picard_config_rule_default.copy()}}

update_config(config_default, config)
config = config_default


rule picard_reorder_sam:
    """Picard: reorder bam file"""
    params: cmd = config['picard']['cmd'] + REORDER_SAM,
            options = config['picard']['reorder_sam']['options'],
            runtime = config['picard']['reorder_sam']['runtime'],
    input: "{prefix}.bam", config['picard']['ref']
    output: "{prefix}.resorted.bam"
    conda: "env.yaml"
    threads: config['picard']['reorder_sam']['threads']
    shell: "{params.cmd} I={input[0]} R={input[1]} O={output} {params.options}"

