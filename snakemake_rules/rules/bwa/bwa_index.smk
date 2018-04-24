# -*- snakemake -*-
include: "bwa.settings"

config_default = {'bwa' :{'bwa_index' : _bwa_config_rule_default.copy()}}

update_config(config_default, config)
config = config_default


rule bwa_index:
    """bwa index a reference"""
    params: runtime = config['bwa']['bwa_index']['runtime']
    input: indexref = config['bwa']['index'] if config['bwa']['index'] else "{prefix}{ext}"
    output: expand("{{prefix}}{{ext,\.[a-z]+}}{bwaext}",\
            bwaext=config['bwa']['index_ext'])
    threads: config['bwa']['bwa_index']['threads']
    conda: "env.yaml"
    shell: "bwa index {input.indexref}"
