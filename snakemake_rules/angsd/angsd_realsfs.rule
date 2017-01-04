# -*- snakemake -*-
include: 'angsd.settings'
include: 'angsd_dosaf.rule'

config_default = {'angsd' :{'realsfs' : _angsd_config_rule_default.copy()}}

update_config(config_default, config)
config = config_default

rule angsd_realsfs:
    """angsd realSFS: estimate the SFS.
    """
    params: options = config['angsd']['realsfs']['options'],
            runtime = config['angsd']['realsfs']['runtime']
    input: idx = "{prefix}.saf.idx"
    output: sfs = "{prefix}.sfs"
    threads: config['angsd']['realsfs']['threads']
    conda: "env.yaml"
    shell: "realSFS -P {threads} {input.idx} > {output.sfs}"

