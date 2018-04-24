# -*- snakemake -*-
include: "samtools.settings.smk"

config_default = {'samtools' :{'gene_region_from_bam' : _samtools_config_rule_default.copy()}}

update_config(config_default, config)
config = config_default

rule samtools_gene_region_from_bam:
    """Extract gene region from bam file"""
    params: options = config['samtools']['gene_region_from_bam']['options'],
            cmd = config['samtools']['cmd'],
            runtime = config['samtools']['gene_region_from_bam']['runtime']
    input: bam = "{prefix}.bam", bed = "{prefix}.region_{region}.{sfx}.bed",
           bai = "{prefix}.bam.bai"
    output: temp("{prefix}.region_{region}.{sfx}.bam")
    threads: config['samtools']['gene_region_from_bam']['threads']
    conda: "env.yaml"
    shell: "samtools view {params.options} -b -L {input.bed} {input.bam} > {output}"
