# -*- snakemake -*-
include: '../ngs.settings'

config_default = {
    'vcftools' : {
        'cmd' : 'vcftools',
        'runtime' : '01:00:00',
    },
}

update_config(config_default, config)
config = config_default


_vcftools_config_rule_default = {
    'options' : '',
    'runtime' : config['vcftools']['runtime'],
    'threads' : 1,
}
