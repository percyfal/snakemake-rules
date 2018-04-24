# -*- snakemake -*-
#
include: "../ngs.settings.smk"

config_default = {
    'ucsc' : { 
        'ref' : config['ngs.settings']['db']['ref'],
        'index' : "",
        'urldownload' : 'http://hgdownload.cse.ucsc.edu/goldenPath/',
        'runtime' : "01:00:00",
    },
}

update_config(config_default, config)
config = config_default

_ucsc_config_rule_default =  {
    'options' : '',
    'runtime' : config['ucsc']['runtime'],
    'threads' : 1
}
