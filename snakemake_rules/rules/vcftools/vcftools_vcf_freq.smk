# -*- snakemake -*-
include: "vcftools.settings.smk"

config_default = {'vcftools' :{'vcf_freq' : _vcftools_config_rule_default.copy()}}

update_config(config_default, config)
config = config_default

rule vcftools_vcf_freq:
    """vcftools: calculate the variant frequencies in a vcf file"""
    params: cmd=config['vcftools']['cmd'],
            options=config['vcftools']['vcf_freq']['options'],
            runtime = config['vcftools']['vcf_freq']['runtime']
    input: "{prefix}.vcf"
    output: "{prefix}.frq"
    threads: config['vcftools']['vcf_freq']['threads']
    shell: "{params.cmd} --freq {params.options} --vcf  {input} --out {wildcards.prefix}"
