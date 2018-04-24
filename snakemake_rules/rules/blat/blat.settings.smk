# -*- snakemake -*-
"""
`blat <https://genome.ucsc.edu/goldenpath/help/blatSpec.html>`_
----------------------------------------------------------------

Rules for blat aligner.

"""

include: '../ngs.settings'

# Programs
BLAT_FATOTWOBIT = 'faToTwoBit'

config_default = {
    'blat' : {
        'ref': config['ngs.settings']['db']['ref'],
        'index': "",
        'runtime' : "01:00:00",
    },
}

update_config(config_default, config)
config = config_default


_blat_config_rule_default = {
    'options' : '',
    'runtime' : config['blat']['runtime'],
    'threads' : 1,
}
