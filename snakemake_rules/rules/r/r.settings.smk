# -*- snakemake -*-
include: '../ngs.settings.smk'

config_default = {
    'r' : {
        'cmd' : 'R',
        'script' : 'Rscript',
        'runtime' : "01:00:00",
        'X11' : True,
    },
}

update_config(config_default, config)
config = config_default


_r_config_rule_default = {
    'options' : '',
    'runtime' : config['r']['runtime'],
    'threads' : 1,
}
