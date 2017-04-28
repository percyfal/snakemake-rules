# -*- snakemake -*-
include: "snpeff.settings"

config_default = {'snpeff' :{'download_database' : _snpeff_config_rule_default.copy()}}

update_config(config_default, config)
config = config_default


rule snpeff_download_database:
    """Download snpEff database. FIXME: fix output so that rule doesn't rerun everytime."""
    params: home=config['snpeff']['home'],
            cmd=config['snpeff']['cmd'],
            runtime = config['snpeff']['download_database']['runtime'],
    output: config['snpeff']['dblist']
    threads: config['snpeff']['download_database']['threads']
    conda: "env.yaml"
    shell: "cd {params.home} && {params.cmd} download {output}"

