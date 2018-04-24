# -*- snakemake -*-
include: '../ngs.settings.smk'

config_default = {
    'bcftools' : {
        'cmd' : 'bcftools',
        'ref' : config['ngs.settings']['db']['ref'],
        'runtime' : '03:00:00',
        'threads' : config['settings']['threads'],
        'options' : {
            'main' : "",
        },
    },
}

update_config(config_default, config)
config = config_default

_bcftools_config_rule_default = {
    'options' : '',
    'runtime' : config['bcftools']['runtime'],
    'threads' : 1,
}
