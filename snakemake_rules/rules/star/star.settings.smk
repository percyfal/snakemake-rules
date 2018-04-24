# -*- snakemake -*-
include: '../ngs.settings'
include: '../comp/comp.settings'

config_default = { 
    'star' : {
        'cmd' : "STAR",
        'ref' : config['ngs.settings']['db']['ref'],
        'extra_ref' : config['ngs.settings']['db']['extra_ref'],
        'index' : os.path.join(os.path.dirname(config['ngs.settings']['db']['ref']), "Genome"),
        'runtime' : "01:00:00",
    },
}

update_config(config_default, config)
config = config_default
