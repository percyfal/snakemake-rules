# -*- snakemake -*-
include: '../main.settings'

config_default = {
    'pandoc' : {
        'cmd' : 'pandoc',
    },
}

update_config(config_default, config)
config = config_default
