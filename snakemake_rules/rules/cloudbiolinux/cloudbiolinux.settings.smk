# -*- snakemake -*-
#
include: "../ngs.settings.smk"

config_default = {
    'cloudbiolinux' : { 
        'ref' : config['ngs.settings']['db']['ref'],
      },
}

update_config(config_default, config)
config = config_default

