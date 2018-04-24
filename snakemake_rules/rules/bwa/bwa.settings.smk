# -*- snakemake -*-
include: '../ngs.settings.smk'

config_default = {
    'bwa' : {
    	'ref' : config['ngs.settings']['db']['ref'],
        'index' : config['ngs.settings']['db']['ref'],
        'index_ext' : ['.amb', '.ann', '.bwt', '.pac', '.sa'],
        'threads' : config['settings']['threads'],
        'runtime' : "01:00:00",
    },
}

update_config(config_default, config)
config = config_default


_bwa_config_rule_default = {
    'options' : '',
    'runtime' : config['bwa']['runtime'],
    'threads' : 1,
}
