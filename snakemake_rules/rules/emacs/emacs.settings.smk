# -*- snakemake -*-
include: '../main.settings.smk'

config_default = {
    'emacs' : {
        'cmd' : 'emacs',
    },
}

update_config(config_default, config)
config = config_default
