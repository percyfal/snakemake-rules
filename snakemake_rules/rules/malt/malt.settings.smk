# -*- snakemake -*-
include: '../ngs.settings'

config_default = {
    'malt': {
        'index': '',
        'ref': config['ngs.settings']['db']['ref'],
        'runtime': "01:00:00",
        'X11': True,
    },
}

update_config(config_default, config)
config = config_default


_malt_config_rule_default = {
    'options' : '',
    'runtime' : config['malt']['runtime'],
    'threads' : 1,
}
