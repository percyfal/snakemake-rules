# -*- snakemake -*-
include: "picard.settings.smk"

config_default = {'picard' :{'calculate_hs_metrics_regions' : _picard_config_rule_default.copy()}}
config_default['picard']['calculate_hs_metrics_regions'].update(
    {
        'targets' : [],
        'bait_regions' : config['ngs.settings']['sequence_capture']['bait_regions'],
        'target_regions' : config['ngs.settings']['sequence_capture']['target_regions'],
    })


update_config(config_default, config)
config = config_default


rule picard_calculate_hs_metrics_regions:
    """Picard: calculate hybrid selection metrics based on regions"""
    params: cmd = config['picard']['cmd'] + CALCULATE_HS_METRICS,
            options = config['picard']['calculate_hs_metrics_regions']['options'],
            runtime = config['picard']['calculate_hs_metrics_regions']['runtime']
    input: bam = "{prefix}.region_{gene}.{sfx}.bam", targets="{prefix}.region_{gene}.targets.dict", baits="{prefix}.region_{gene}.baits.dict", ref=config['picard']['ref']
    output: metrics = "{prefix}.region_{gene}.{sfx}.hs_metrics"
    conda: "env.yaml"
    threads: config['picard']['calculate_hs_metrics_regions']['threads']
    shell: "{params.cmd} {params.options} TI={input.targets} BI={input.baits} I={input.bam} O={output.metrics} R={input.ref}"

