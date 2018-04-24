# -*- snakemake -*-
include: '../ngs.settings.smk'

config_default = { 
    'sga' : {
        'cmd' : 'sga',
        'threads' : config['settings']['threads'],
        'options' : "",
        'runtime' : "01:00:00",
    },
}

update_config(config_default, config)
config = config_default

_sga_config_rule_default = {
    'options' : '',
    'runtime' : config['sga']['runtime'],
    'threads' : 1,
}

