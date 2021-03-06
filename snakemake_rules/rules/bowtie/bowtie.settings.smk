# -*- snakemake -*-
"""
`bowtie <http://bowtie-bio.sourceforge.net/index.shtml>`_
------------------------------------------------------------

Bowtie - An ultrafast memory-efficient short read aligner

"""
include: '../ngs.settings.smk'
include: '../samtools/samtools.settings.smk'

config_default = {
    'bowtie' : {
        'build_ext': [".1.ebwt", ".2.ebwt", ".3.ebwt",
                      ".4.ebwt", ".rev.1.ebwt", ".rev.2.ebwt"],
        'cmd': 'bowtie',
        'index' : os.path.splitext(config['ngs.settings']['db']['ref'])[0],
        'ref' : config['ngs.settings']['db']['ref'],
        'rg_fn' : None,
        'runtime' : "01:00:00",
        'threads' : config['settings']['threads'],
    },
}

update_config(config_default, config)
config = config_default

_bowtie_config_rule_default = {
    'options' : '',
    'runtime' : config['bowtie']['runtime'],
    'threads' : 1,
}
