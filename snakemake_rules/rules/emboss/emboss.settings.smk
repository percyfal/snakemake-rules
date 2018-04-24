# -*- snakemake -*-
include: '../ngs.settings.smk'


config_default = { 
    'emboss' : {
        'runtime' : "01:00:00",
    },
}

update_config(config_default, config)
config = config_default


_emboss_config_rule_default = {
    'options' : '',
    'runtime' : config['emboss']['runtime'],
    'threads' : 1,
}
