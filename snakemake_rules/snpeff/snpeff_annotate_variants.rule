# -*- snakemake -*-
include: "snpeff.settings"

config_default = {'snpeff' :{'annotate' : _snpeff_config_rule_default.copy()}}

update_config(config_default, config)
config = config_default

cmd = re.sub("-Xmx[0-9a-zA-Z]+", "-Xmx{mem}".format(mem=config['snpeff']['annotate']['java_mem']), config['snpeff']['cmd'])


rule snpeff_annotate_variants:
    """snpEff: annotate variants"""
    params: cmd = " ".join([config['snpeff']['cmd'], "ann"]),
            options=config['snpeff']['annotate']['options'],
            genome_version=config['snpeff']['genome_version'],
            runtime=config['snpeff']['annotate']['runtime'],
    input: "{prefix}.vcf"
    output: "{prefix}.annotated.{sfx}"
    threads: config['snpeff']['annotate']['threads']
    conda: "env.yaml"
    shell: "{params.cmd} ann {params.options} {params.genome_version} {input} -i vcf -o {wildcards.sfx} > {output}"
