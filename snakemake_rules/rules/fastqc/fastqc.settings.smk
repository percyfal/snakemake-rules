# -*- snakemake -*-
include: "../ngs.settings.smk"

config_default = {
    'fastqc' : {
        'cmd' : "fastqc",
        'options' : "",
        'runtime' : "00:20:00",
        'threads' : config['settings']['threads'],
    },
}

update_config(config_default, config)
config = config_default
