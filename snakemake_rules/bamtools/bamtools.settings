# -*- snakemake -*-
include: '../ngs.settings'

config_default = { 
    'bamtools' : {
        'ref' : config['ngs.settings']['db']['ref'],
        'cmd' : "bamtools",
        'home' : "",
        'threads' : config['settings']['threads'],
        'options' : "",
        'runtime' : "01:00:00",
    },
}

update_config(config_default, config)
config = config_default


_bamtools_config_rule_default = {
    'options' : {'mapQuality' : ">=255"},
    'regions' : config['ngs.settings']['regions'],
    'runtime' : config['bamtools']['runtime'],
    'threads' : 1,
}
