# -*- snakemake -*-
include: '../ngs.settings'

config_default = {
    'freebayes' : {
        'cmd' : 'freebayes',
        'ref' : config['ngs.settings']['db']['ref'],
        'runtime' : '100:00:00',
        'options' : "",
        'threads' : 1,
        'targets' : config['ngs.settings']['regions'],
    },
}

update_config(config_default, config)
config = config_default

_freebayes_config_rule_default = {
    'ref': config['freebayes']['ref'],
    'runtime': config['freebayes']['runtime'],
    'threads': 1,
}
