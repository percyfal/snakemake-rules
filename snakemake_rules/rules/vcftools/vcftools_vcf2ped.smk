# -*- snakemake -*-
include: "vcftools.settings.smk"

config_default = {'vcftools' :{'vcf2ped' : _vcftools_config_rule_default.copy()}}

update_config(config_default, config)
config = config_default

rule vcftools_vcf2ped:
    """vcftools: convert ped file to vcf."""
    params: cmd=config['vcftools']['cmd'],
            options=config['vcftools']['vcf2ped']['options'],
            runtime = config['vcftools']['vcf2ped']['runtime']
    input: vcf = "{prefix}.vcf"
    output: ped = "{prefix}.ped"
    threads: config['vcftools']['vcf2ped']['threads']
    shell: "out={output.ped} && {params.cmd} {params.options} --vcf {input.vcf} --plink --out ${{out%.ped}}"

