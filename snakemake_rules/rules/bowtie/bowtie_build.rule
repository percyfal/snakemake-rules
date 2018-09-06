# -*- snakemake -*-
include: "bowtie.settings"

config_default = {'bowtie' :{'bowtie_build' : _bowtie_config_rule_default.copy()}}
config_default['bowtie']['bowtie_build'].update({'cmd': 'bowtie-build'})

update_config(config_default, config)
config = config_default

rule bowtie_build:
    """Bowtie build index"""
    params: cmd = config['bowtie']['bowtie_build']['cmd'],
            runtime = config['bowtie']['bowtie_build']['runtime'],
            options = config['bowtie']['bowtie_build']['options'],
    input: config['bowtie']['ref'] if config['bowtie']['ref'] else "{prefix}"
    output: expand("{{prefix}}{ext}", ext=config['bowtie']['build_ext'])
    threads: config['bowtie']['bowtie_build']['threads']
    conda: "env.yaml"
    shell: "{params.cmd} {params.options} --threads {threads} {input} {wildcards.prefix}"
