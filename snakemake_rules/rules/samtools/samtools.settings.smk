# -*- snakemake -*-
include: '../ngs.settings.smk'

config_default = { 
    'samtools' : {
        'cmd' : 'samtools',
        'ref' : config['ngs.settings']['db']['ref'],
        'threads' : config['settings']['threads'],
        'options' : "",
        'runtime' : "01:00:00",
    },
}

update_config(config_default, config)
config = config_default

_samtools_config_rule_default = {
    'options' : '',
    'runtime' : config['samtools']['runtime'],
    'threads' : 1,
}

