# -*- snakemake -*-
# Qualimap: http://qualimap.bioinfo.cipf.es/
#
include: '../ngs.settings.smk'

config_default = {
    'qualimap' : {
        'cmd' : 'qualimap',
        'threads' : config['settings']['threads'],
        'java_mem' : config['settings']['java']['java_mem'],
    },
}


update_config(config_default, config)
config = config_default

