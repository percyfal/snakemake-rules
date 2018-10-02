# -*- snakemake -*-
include: "picard.settings.smk"
include: "picard_build_bam_index.smk"

config_default = {'picard' :{'calculate_hs_metrics' : _picard_config_rule_default.copy()}}
config_default['picard']['calculate_hs_metrics'].update(
    {
        'targets' : [],
        'bait_regions' : config['ngs.settings']['sequence_capture']['bait_regions'],
        'target_regions' : config['ngs.settings']['sequence_capture']['target_regions'],
    })


update_config(config_default, config)
config = config_default

rule picard_calculate_hs_metrics:
    """Picard: calculate hybrid selection metrics"""
    params: cmd = config['picard']['cmd'] + CALCULATE_HS_METRICS,
            options = config['picard']['calculate_hs_metrics']['options'],
            suffix = ".hs_metrics",
            runtime = config['picard']['calculate_hs_metrics']['runtime'],
    input: bam = "{prefix}.bam",
           bai = "{prefix}.bai",
           ref = config['picard']['ref'],
           bait_regions = config['picard']['calculate_hs_metrics']['bait_regions'],
           target_regions = config['picard']['calculate_hs_metrics']['target_regions']
    output: metrics = "{prefix}.hs_metrics"
    conda: "env.yaml"
    threads: config['picard']['calculate_hs_metrics']['threads']
    shell: "{params.cmd} {params.options} TI={input.target_regions} BI={input.bait_regions} I={input.bam} O={output.metrics} R={input.ref}"

