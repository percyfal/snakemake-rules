# -*- snakemake -*-
include: 'gatk.settings'
include: 'gatk_base_recalibrator.rule'

config_default = {'gatk' : {'print_reads' : _gatk_config_rule_default.copy()}}

update_config(config_default, config)
config = config_default


rule gatk_print_reads:
    """Run GATK PrintReads"""
    params: cmd = config['gatk']['cmd'] + " -T " + PRINT_READS,
            options = config['gatk']['print_reads']['options'],
            runtime = config['gatk']['print_reads']['runtime'],
    input: bam = "{prefix}.bam",
           bai = "{prefix}.bai",
           recal = "{prefix}.recal_data.grp",
           ref = config['gatk']['print_reads']['ref']
    output: bam = "{prefix}.recal.bam", bai = "{prefix}.recal.bai"
    threads: config['gatk']['print_reads']['threads']
    conda: "env.yaml"
    shell: "{params.cmd} {params.options} -R {input.ref} -I {input.bam} -BQSR {input.recal} -o {output.bam}"
