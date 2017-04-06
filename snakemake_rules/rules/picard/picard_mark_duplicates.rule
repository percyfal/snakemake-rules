# -*- snakemake -*-
include: "picard.settings"

config_default = {'picard' :{'mark_duplicates' : _picard_config_rule_default.copy()}}
config_default['picard']['mark_duplicates'].update({'targets' : []})

update_config(config_default, config)
config = config_default


rule picard_mark_duplicates:
    """Picard: mark duplicates"""
    params: cmd = config['picard']['cmd'] + MARK_DUPLICATES,
            options = config['picard']['mark_duplicates']['options'],
            runtime = config['picard']['mark_duplicates']['runtime']
    input: bam = "{prefix}.bam", bai = "{prefix}.bai"
    output: bam = "{prefix}.dup.bam", metrics = "{prefix}.dup.dup_metrics"
    threads: config['picard']['mark_duplicates']['threads']
    shell: "{params.cmd} I={input.bam} O={output.bam} {params.options} M={output.metrics}"



rule picard_mark_duplicates_log:
    """Picard: mark duplicates, log output"""
    params: runtime = "00:10:00"
    input: metrics = "{prefix}.dup.dup_metrics"
    output: metrics = "{prefix}.dup.dup_metrics.log"
    conda: "env.yaml"
    threads: 1
    shell: "cp {input.metrics} {output.metrics}"

localrules: picard_mark_duplicates_log           
