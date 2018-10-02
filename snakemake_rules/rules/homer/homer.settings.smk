# -*- snakemake -*-
#
# HOMER: Software for motif discovery and next-gen sequencing analysis
#
# http://homer.salk.edu/homer/motif/
#
include: '../ngs.settings.smk'

config_default = {
    'homer' : {
        'options' : '',
        'genome_build': 'hg38',
        'runtime' : '01:00:00',
    },
}

update_config(config_default, config)
config = config_default


_homer_config_rule_default = {
    'options' : '',
    'runtime' : config['homer']['runtime'],
    'threads' : 1,
}
