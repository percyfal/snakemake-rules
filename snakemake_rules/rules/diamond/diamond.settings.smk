# -*- snakemake -*-
include: '../ngs.settings'

config_default = {
    'diamond' : {
        'ref' : config['ngs.settings']['db']['ref'],
        'options' : "",
    },
}

update_config(config_default, config)
config = config_default
