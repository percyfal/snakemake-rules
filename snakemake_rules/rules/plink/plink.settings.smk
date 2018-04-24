# -*- snakemake -*-
include: "../ngs.settings.smk"

config_default = { 
    'plink' : {
        'common_options' : "--noweb",
        'cmd' : "plink2",
        'chr' : "1",
        'runtime' : "01:00:00",
        'threads' : 1,
    },
}

update_config(config_default, config)
config = config_default

config_default2 = {
    'plink' : {
        'options' : config['plink']['common_options']
    },
}

update_config(config_default2, config)
config = config_default2
