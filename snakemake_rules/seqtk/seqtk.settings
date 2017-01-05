# -*- snakemake -*-
include: '../ngs.settings'


config_default = { 
    'seqtk' : {
        'cmd' : "seqtk",
        'runtime' : '01:00:00',
        'options' : "",
    },
}

update_config(config_default, config)
config = config_default


_seqtk_config_rule_default = {
    'options' : '',
    'runtime' : config['seqtk']['runtime'],
    'threads' : 1,
}
