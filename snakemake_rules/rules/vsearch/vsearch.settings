# -*- snakemake -*-
include: "../ngs.settings"

config_default = {
    'vsearch' : {
        'cmd' : "vsearch",
        'options' : "",
        'threads' : config['settings']['threads'],
    },
}

update_config(config_default, config)
config = config_default
