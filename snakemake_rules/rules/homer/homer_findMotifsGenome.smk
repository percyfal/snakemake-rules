# -*- snakemake -*-
include: "homer.settings.smk"

config_default = {'homer' :{'findMotifsGenome' : _homer_config_rule_default.copy()}}
config_default['homer']['findMotifsGenome'].update({'size': 10})

update_config(config_default, config)
config = config_default

rule findMotifsGenome:
    """homer: run findMotifsGenome.pl for finding enriched motifs in
    genomic regions"""
    params: cmd = "findMotifsGenome.pl",
            build = config['homer']['genome_build'],
            size = config['homer']['findMotifsGenome']['size'],
            options = " ".join([config['homer']['findMotifsGenome']['options']]),
            runtime = config['homer']['findMotifsGenome']['runtime'],
    input: bed = "{prefix}.bed"
    output: homermotifs = "{prefix}.findMotifsGenome/homerMotifs.all.motifs"
    log: "{prefix}.findMotifsGenome.log"
    threads: config['homer']['findMotifsGenome']['threads']
    shell: "{params.cmd} {input.bed} {params.build} $(dirname {output.homermotifs}) -size {params.size} {params.options} 2> {log}"
