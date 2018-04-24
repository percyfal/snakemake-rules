# -*- snakemake -*-
include: "ucsc.settings.smk"

config_default = {'ucsc': {'genePredToBed' : _ucsc_config_rule_default.copy()}}
config_default['ucsc']['genePredToBed'].update({'cmd' : 'genePredToBed'})

update_config(config_default, config)
config = config_default


rule ucsc_genePredToBed:
    """Run genePredToBed"""
    params: cmd = config['ucsc']['genePredToBed']['cmd'],
            options = config['ucsc']['genePredToBed']['options'],
            runtime = config['ucsc']['genePredToBed']['runtime']
    input: genepred = "{prefix}.genePred"
    output: bed = "{prefix}.bed"
    threads: config['ucsc']['genePredToBed']['threads']
    conda: "envs/ucsc_genepredtobed.yaml"
    shell: "{params.cmd} {params.options} {input.genepred} {output.bed}"
