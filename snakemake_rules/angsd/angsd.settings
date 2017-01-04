# -*- snakemake -*-
include: '../ngs.settings'

config_default = {
    'angsd' : {
        'cmd' : 'angsd',
        'ref' : config['ngs.settings']['db']['ref'],
        'anc' : config['ngs.settings']['db']['ref'],
        'threads' : config['settings']['threads'],
        'options' : '',
        'gl' : 1,
        'majorminor' : 1,
        'nChr' : 1, # Number of chromosomes
        'runtime' : '01:00:00',
    },
}

update_config(config_default, config)
config = config_default


_angsd_config_rule_default = {
    'options' : '',
    'runtime' : config['angsd']['runtime'],
    'threads' : 1,
}
