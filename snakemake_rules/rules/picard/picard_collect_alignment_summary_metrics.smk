# -*- snakemake -*-
include: "picard.settings.smk"
include: "picard_build_bam_index.smk"

config_default = {'picard' :{'collect_alignment_summary_metrics' : _picard_config_rule_default.copy()}}
config_default['picard']['collect_alignment_summary_metrics'].update({'targets' : []})

update_config(config_default, config)
config = config_default


rule picard_collect_alignment_summary_metrics:
    """Picard: collect alignment summary metrics"""
    params: cmd = config['picard']['cmd'] + COLLECT_ALIGNMENT_SUMMARY_METRICS,
            options = config['picard']['collect_alignment_summary_metrics']['options'],
            suffix = ".align_metrics",
            runtime = config['picard']['collect_alignment_summary_metrics']['runtime']
    input: bam = "{prefix}.bam", bai = "{prefix}.bai", ref = config['picard']['ref']
    output: metrics = "{prefix}.align_metrics"
    conda: "env.yaml"
    threads: config['picard']['collect_alignment_summary_metrics']['threads']
    shell: "{params.cmd} {params.options} I={input.bam} O={output.metrics} R={input.ref}"

