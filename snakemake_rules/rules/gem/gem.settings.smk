# -*- snakemake -*-
"""
`GEM <http://algorithms.cnag.cat/wiki/The_GEM_library>`_
----------------------------------------------------------

The GEM library
"""
include: "../ngs.settings.smk"

config_default = {
    'gem' : {
        'runtime': "01:00:00",
    },
}

update_config(config_default, config)
config = config_default


_gem_config_rule_default = {
    'options' : '',
    'runtime' : config['gem']['runtime'],
    'threads' : 1,
}
