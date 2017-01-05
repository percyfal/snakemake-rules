# -*- snakemake -*-
#
# Install with pip; use recent version of gcc (>=4.9)
#
include: "../ngs.settings"

config_default = {
    'macs2' : {
        'cmd' : 'macs2',
        'runtime' : "01:00:00",
    },
}

update_config(config_default, config)
config = config_default


_macs2_config_rule_default = {
    'options' : '',
    'runtime' : config['macs2']['runtime'],
    'threads' : 1,
}
