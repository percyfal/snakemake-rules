# -*- snakemake -*-
include: "../main.settings"


config_default = {
    'jellyfish' : {
        'cmd': 'jellyfish',
        'threads': config['settings']['threads'],
        'runtime': "01:00:00",
    },
}

update_config(config_default, config)
config = config_default
