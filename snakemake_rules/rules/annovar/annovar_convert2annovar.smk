# -*- snakemake -*-
include: "annovar.settings.smk"

config_default = {'annovar' :{'convert2annovar' : _annovar_config_rule_default.copy()}}
config_default['annovar']['convert2annovar'].update(
    {
        'options' : "--includeinfo",
        'vcf_format' : "vcf4",
    })

update_config(config_default, config)
config = config_default


rule annovar_convert2annovar:
    """Convert data to annovar format"""
    params: cmd='convert2annovar.pl',
            options=config['annovar']['convert2annovar']['options'],
            vcf_format = config['annovar']['convert2annovar']['vcf_format'],
            runtime = config['annovar']['convert2annovar']['runtime']
    input: "{prefix}.vcf"
    output: "{prefix}.avinput"
    threads: config['annovar']['convert2annovar']['threads']
    shell: "{params.cmd} {params.options} -format {params.vcf_format} {input} > {output}"
