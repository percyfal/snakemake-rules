# -*- snakemake -*-
include: '../ngs.settings'

config_default = {
    'htslib' : {
        'threads' : config['settings']['threads'],
        'runtime': "01:00:00",
    },
}

update_config(config_default, config)
config = config_default

_htslib_config_rule_default = {
    'options' : '',
    'runtime' : config['htslib']['runtime'],
    'threads' : 1,
}
