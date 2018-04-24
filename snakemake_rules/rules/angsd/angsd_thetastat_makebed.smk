# -*- snakemake -*-
include: 'angsd.settings.smk'
include: 'angsd_dothetas.smk'

config_default = {'angsd' :{'thetastat_makebed' : _angsd_config_rule_default.copy()}}

update_config(config_default, config)
config = config_default

rule angsd_thetastat_makebed:
    """angsd thetastat_makebed: perform neutrality test statistic on thetas.gz output

    makebed: generate input for dostat
    """
    params: options = config['angsd']['thetastat_makebed']['options'],
            runtime = config['angsd']['thetastat_makebed']['runtime']
    threads: config['angsd']['thetastat_makebed']['threads']
    input: thetasgz = "{prefix}.thetas.gz"
    output: bin = "{prefix}.bin", idx = "{prefix}.idx"
    conda: "env.yaml"
    shell: "thetaStat make_bed {input.thetasgz} {wildcards.prefix}"
