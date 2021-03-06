# -*- snakemake -*-
include: "bowtie2.settings.smk"

config_default = {'bowtie2' :{'bowtie2_build' : _bowtie2_config_rule_default.copy()}}
config_default['bowtie2']['bowtie2_build'].update({'cmd': 'bowtie2-build'})

update_config(config_default, config)
config = config_default

rule bowtie2_build:
    """Bowtie build index"""
    params: cmd = config['bowtie2']['bowtie2_build']['cmd'],
            runtime = config['bowtie2']['bowtie2_build']['runtime'],
            options = config['bowtie2']['bowtie2_build']['options'],
    input: config['bowtie2']['ref'] if config['bowtie2']['ref'] else "{prefix}"
    output: expand("{{prefix}}{ext}", ext=config['bowtie2']['build_ext'])
    threads: config['bowtie2']['bowtie2_build']['threads']
    shell: "{params.cmd} {params.options} --threads {threads} {input} {wildcards.prefix}"
