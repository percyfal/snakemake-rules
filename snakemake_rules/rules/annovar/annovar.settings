# -*- snakemake -*-
include: '../ngs.settings'

config_default = {
    'annovar' : {
        'buildver' : config['ngs.settings']['db']['ref'],
        'dblist' : ["dgv", "genomicSuperDups", "gwascatalog", "tfbs", "wgEncodeRegTfbsClustered", "wgEncodeRegDnaseClustered", "phastConsElements46way"],
        'dblist_webfrom_annovar' : ["1000g2012apr", "cosmic64", "esp6500si_all", "esp6500si_ea", "ljb_all", "snp137", "refGene", "avsift"],
        'home': os.getenv("ANNOVAR_HOME", os.curdir),
        'runtime' : '01:00:00',
    },
}

update_config(config_default, config)
config = config_default

config_default2 = {
    'annovar': {
        'db' : os.path.join(config['annovar']['home'], "/humandb")
    }
}

update_config(config_default2, config)
config = config_default2


_annovar_config_rule_default = {
    'options' : '',
    'runtime' : config['annovar']['runtime'],
    'threads' : 1,
}
