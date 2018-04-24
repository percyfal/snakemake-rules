# -*- snakemake -*-
include: '../ngs.settings.smk'

config_default = {
    'vcflib' : {
        'ref' : config['ngs.settings']['db']['ref'],
        'runtime' : '01:00:00',
    },
}

update_config(config_default, config)
config = config_default

_vcflib_config_rule_default = {
    'options' : '',
    'runtime' : config['vcflib']['runtime'],
    'threads' : 1,
}
