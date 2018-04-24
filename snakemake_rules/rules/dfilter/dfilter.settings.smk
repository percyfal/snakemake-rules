# -*- snakemake -*-
#
# http://collaborations.gis.a-star.edu.sg/~cmb6/kumarv1/dfilter/
#
# Requires manual install
#
include: '../ngs.settings.smk'


config_default = {
    'dfilter' : {
        'ref' : config['ngs.settings']['db']['ref'],
        'options' : "-lpval=2 -ks=50 -bs=100 -wig",
        'cmd' : 'run_dfilter.sh',
    },
}

update_config(config_default, config)
config = config_default
