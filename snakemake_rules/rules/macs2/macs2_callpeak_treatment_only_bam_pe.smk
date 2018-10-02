# -*- snakemake -*-
include: "macs2.settings.smk"

config_default = {'macs2' :{'callpeak' : _macs2_config_rule_default.copy()}}
config_default['macs2']['callpeak'].update({'options': '-B'})


update_config(config_default, config)
config = config_default
 

rule macs2_callpeak_treatment_only_bam_pe:
    """MACS callpeak on sample with no control, paired end."""
    params: cmd = config['macs2']['cmd'],
            options = config['macs2']['callpeak']['options'],
            runtime = config['macs2']['callpeak']['runtime']
    input: bam = "{prefix}.bam"
    output: xls = "{prefix}_peaks.xls",
            bedGraph = ["{prefix}_control_lambda.bdg", "{prefix}_treat_pileup.bdg"] if "-B" in config['macs2']['callpeak']['options'] else [],
            summits = "{prefix}_summits.bed",
            narrowPeak = "{prefix}_peaks.narrowPeak",
    log: "{prefix}.macs2.log"
    threads: config['macs2']['callpeak']['threads']
    conda: "env.yaml"
    shell: "{params.cmd} callpeak {params.options} -n {wildcards.prefix} -t {input.bam} -f BAMPE &> {log}"
