# -*- snakemake -*-
include: "ucsc.settings"

config_default = {'ucsc': {'gtfToGenePred' : _ucsc_config_rule_default.copy()}}
config_default['ucsc']['gtfToGenePred'].update(
    {
        'cmd' : 'gtfToGenePred',
        'options' : '-genePredExt -ignoreGroupsWithoutExons',
     })


update_config(config_default, config)
config = config_default


rule ucsc_gtfToGenepred:
    """Run gtfToGenePred"""
    params: cmd = config['ucsc']['gtfToGenePred']['cmd'],
            options = config['ucsc']['gtfToGenePred']['options'],
            runtime = config['ucsc']['gtfToGenePred']['runtime']
    input: gtf = "{prefix}.gtf"
    output: genepred = "{prefix}.genePred"
    threads: config['ucsc']['gtfToGenePred']['threads']
    conda:
        join("envs", "ucsc_gtftogenepred.yaml")
    shell: "{params.cmd} {params.options} {input.gtf} {output.genepred}"
