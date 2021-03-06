# -*- snakemake -*-
include: 'gatk.settings.smk'
include: '../picard/picard_build_bam_index.smk'

config_default = {'gatk': {'clip_reads' : _gatk_config_rule_default.copy()}}

update_config(config_default, config)
config = config_default


rule gatk_clip_reads:
    """Run GATK ClipReads"""
    params: cmd = config['gatk']['cmd'] + " -T " + CLIP_READS,
            options = " ".join([config['gatk']['clip_reads']['options']]),
            runtime = config['gatk']['clip_reads']['runtime']
    input: bam = "{prefix}.bam", bai = "{prefix}.bai", ref = config['gatk']['ref']
    output: "{prefix}.clip.bam"
    threads: config['gatk']['clip_reads']['threads']
    conda: "env.yaml"
    shell: "{params.cmd} {params.options} -R {input.ref} -I {input.bam} -o {output}"

