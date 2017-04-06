# -*- snakemake -*-
include: "picard.settings"

config_default = {'picard' :{'build_bam_index' : _picard_config_rule_default.copy()}}

update_config(config_default, config)
config = config_default


rule picard_build_bam_index:
    """Picard: build bam index from bam file"""
    params: cmd = config['picard']['cmd'] +  BUILD_BAM_INDEX,
            options = config['picard']['build_bam_index']['options'],
            runtime = config['picard']['build_bam_index']['runtime']
    input: "{prefix}.bam"
    output: "{prefix}.bai"
    conda: "env.yaml"
    threads: config['picard']['build_bam_index']['threads']
    shell: "{params.cmd} I={input} O={output} {params.options}"


