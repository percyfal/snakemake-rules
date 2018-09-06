# -*- snakemake -*-
include: "ucsc.settings.smk"

rule ucsc_pseudo:
    """Pseudo rule; eliminates circular rule dependency of autosome.fa -> fa -> autosome.fa"""
    output: ref = protected(config['ucsc']['ref'])

localrules: ucsc_pseudo            
