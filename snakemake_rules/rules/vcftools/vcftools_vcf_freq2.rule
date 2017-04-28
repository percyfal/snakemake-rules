# -*- snakemake -*-
include: "vcftools.settings"

config_default = {'vcftools' :{'vcf_freq2' : _vcftools_config_rule_default.copy()}}

update_config(config_default, config)
config = config_default


rule vcftools_vcf_freq2:
    """vcftools: calculate the variant frequencies in a vcf file using option --freq2"""
    params: cmd=config['vcftools']['cmd'],
            options=config['vcftools']['vcf_freq2']['options'],
            runtime = config['vcftools']['vcf_freq2']['runtime']
    input: "{prefix}.vcf"
    output: "{prefix}.2.frq"
    threads: config['vcftools']['vcf_freq2']['threads']
    shell: "{params.cmd} --freq2 {params.options} --vcf {input} --out {wildcards.prefix}.2"

