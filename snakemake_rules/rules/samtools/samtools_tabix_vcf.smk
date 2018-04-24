# -*- snakemake -*-
include: "samtools.settings.smk"

config_default = {'samtools' :{'tabix_vcf' : _samtools_config_rule_default.copy()}}
config_default['samtools']['tabix_vcf'].update({'cmd' : 'tabix'})

update_config(config_default, config)
config = config_default


rule samtools_tabix_vcf:
    """Run samtools tabix on gzipped vcf file"""
    params: options = "-p vcf " + config['samtools']['tabix_vcf']['options'],
            cmd = config['samtools']['tabix_vcf']['cmd'],
            runtime = config['samtools']['tabix_vcf']['runtime']
    input: "{prefix}.vcf.gz"
    output: "{prefix}.vcf.gz.tbi"
    threads: config['samtools']['tabix_vcf']['threads']
    conda: "env.yaml"
    shell: "{params.cmd} {params.options} {input}"
