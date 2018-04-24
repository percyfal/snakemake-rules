# -*- snakemake -*-
include: "samtools.settings.smk"

config_default = {'samtools' :{'sam2bam' : _samtools_config_rule_default.copy()}}

update_config(config_default, config)
config = config_default


rule samtools_sam2bam:
    """Convert sam file to bam."""
    params: options = config['samtools']['sam2bam']['options'],
            cmd = config['samtools']['cmd'],
            runtime = config['samtools']['sam2bam']['runtime'],
    input: "{prefix}.sam"
    output: "{prefix}.bam"
    threads: config['samtools']['sam2bam']['threads']
    conda: "env.yaml"
    shell: "{params.cmd} view {params.options} -Sb {input} > {output}"

