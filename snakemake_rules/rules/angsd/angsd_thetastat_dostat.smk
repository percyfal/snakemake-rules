# -*- snakemake -*-
include: 'angsd.settings.smk'
include: 'angsd_dothetas.smk'
include: 'angsd_thetastat_makebed.smk'

config_default = {'angsd' :{'thetastat_dostat' : _angsd_config_rule_default.copy()}}

update_config(config_default, config)
config = config_default

rule angsd_thetastat_dostat:
    """angsd thetastat dostat: estimate for every chromosome
    """
    params:
        options = " ".join(['-nChr {}'.format(config['angsd']['nChr']),
                            config['angsd']['thetastat_dostat']['options']]),
        runtime = config['angsd']['thetastat_dostat']['runtime']
    input: bin = "{prefix}.bin", idx = "{prefix}.idx"
    output: thetas = "{prefix}.pestPG"
    threads: config['angsd']['thetastat_dostat']['threads']
    conda: "env.yaml"
    shell: "thetaStat do_stat {wildcards.prefix} {params.options}"
