# -*- snakemake -*-
include: "ucsc.settings.smk"

config_default = {'ucsc': {'download_2bit': _ucsc_config_rule_default.copy()}}

update_config(config_default, config)
config = config_default

rule ucsc_download_2bit:
    """Download 2bit file from ucsc"""
    params: urlinput = os.path.join(config['ucsc']['urldownload'], "{build}", "bigZips", "{build}.2bit"),
            runtime = config['ucsc']['runtime']
    output: os.path.join("{path}", "{build}", "ucsc", "{build}.2bit")
    threads: config['ucsc']['download_2bit']['threads']
    shell: "mkdir -p `dirname {output}` && wget {params.urlinput} -O {output}"

