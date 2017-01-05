# -*- snakemake -*-
include: "ercc.settings"

config_default = {'ercc' :{'download_metadata' : _ercc_config_rule_default.copy()}}

update_config(config_default, config)
config = config_default

rule ercc_download_metadata:
    params: source = config['ercc']['source'],
            runtime = config['ercc']['download_metadata']['runtime']
    output: temp('cms_095047.txt')
    threads: config['ercc']['download_metadata']['threads']
    shell: "wget {params.source} -O {output}"

