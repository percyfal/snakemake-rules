# -*- snakemake -*-
include: '../ngs.settings.smk'

config_default = {
    'mapdamage2' : {
        'ref' : config['ngs.settings']['db']['ref'],
        'runtime' : "01:00:00",
        'threads' : 1,
        'options' : "",
    },
}

update_config(config_default, config)
config = config_default
