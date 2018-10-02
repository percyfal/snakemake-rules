# -*- snakemake -*-
include: "../ngs.settings.smk"

config_default = {
    'multiqc' : {
        'cmd' : "multiqc",
        'options' : "-f",
        'runtime' : '00:60:00',
        'inputs' : [],
        'threads' : 1,
    },
}

update_config(config_default, config)
config = config_default
