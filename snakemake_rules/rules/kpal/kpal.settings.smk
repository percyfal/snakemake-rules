# -*- snakemake -*-
include: '../ngs.settings.smk'

config_default = { 
    'kpal' : {
        'cmd' : "kpal",
        'runtime' : '24:00:00',
        'options' : "",
    },
}

update_config(config_default, config)
config = config_default

_kpal_config_rule_default = {
    'options' : '',
    'runtime' : config['kpal']['runtime'],
    'threads' : 1,
}
