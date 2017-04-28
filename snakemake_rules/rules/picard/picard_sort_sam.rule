# -*- snakemake -*-
include: "picard.settings"

config_default = {'picard' :{'sort_sam' : _picard_config_rule_default.copy()}}
config_default['picard']['sort_sam'].update({'options' : "SORT_ORDER=coordinate"})

update_config(config_default, config)
config = config_default


rule picard_sort_sam:
    """Picard: sort bam file"""
    params: cmd = config['picard']['cmd'] + SORT_SAM,
            options = config['picard']['options'],
            sortsam_options = config['picard']['sort_sam']['options'],
            runtime = config['picard']['sort_sam']['runtime']
    input: "{prefix}.bam"
    output: "{prefix}.sort.bam"
    conda: "env.yaml"
    threads: config['picard']['sort_sam']['threads']
    shell: "{params.cmd} I={input} O={output} {params.options} {params.sortsam_options}"

