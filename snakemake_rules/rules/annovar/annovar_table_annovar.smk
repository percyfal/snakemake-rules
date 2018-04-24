# -*- snakemake -*-
include: "annovar.settings.smk"

config_default = {'annovar' :{'table_annovar' : _annovar_config_rule_default.copy()}}
config_default['annovar']['table_annovar'].update({'options' : "--otherinfo"})

update_config(config_default, config)
config = config_default


rule annovar_table_annovar:
    """annovar: run table_annovar. Currently only defined for hg19."""
    params: cmd='table_annovar.pl',
            options=config['annovar']['table_annovar']['options'],
            db=config['annovar']['db'],
            buildver=config['annovar']['buildver'],
            runtime=config['annovar']['table_annovar']['runtime'],
    input: "{prefix}.avinput"
    output: "{prefix}.avinput.hg19_multianno.txt"
    threads: config['annovar']['table_annovar']['threads']
    shell: "{params.cmd} {input} {params.db} {params.options} --buildver hg19"

