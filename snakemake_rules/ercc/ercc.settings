# -*- snakemake -*-
include: "../ngs.settings"

config_default = {
    'ercc' : {
        'runtime': "01:00:00",
        # Using Lifetech's source for now
        'source': 'https://tools.lifetechnologies.com/content/sfs/manuals/cms_095047.txt',
    },
}

update_config(config_default, config)
config = config_default


_ercc_config_rule_default = {
    'options' : '',
    'runtime' : config['ercc']['runtime'],
    'threads' : 1,
}
