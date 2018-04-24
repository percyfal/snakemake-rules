# -*- snakemake -*-
"""
`bowtie2 <http://bowtie-bio.sourceforge.net/bowtie2/index.shtml>`_
--------------------------------------------------------------------

Bowtie 2 - Fast and sensitive read alignment
"""
include: '../ngs.settings.smk'
include: '../samtools/samtools.settings.smk'


config_default = {
    'bowtie2' : {
        'build_ext': [".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2"],
        'cmd' : 'bowtie2',
        'index' : os.path.splitext(config['ngs.settings']['db']['ref'])[0],
        'ref' : config['ngs.settings']['db']['ref'],
        'rg_fn' : None,
        'runtime' : "01:00:00",
    },
}

update_config(config_default, config)
config = config_default


_bowtie2_config_rule_default = {
    'options' : '',
    'runtime' : config['bowtie2']['runtime'],
    'threads' : 1,
}
