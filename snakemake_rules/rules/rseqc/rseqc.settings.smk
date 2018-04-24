# -*- snakemake -*-
include: '../ngs.settings'
# include: '../misc/conversion.rules'

# Programs
CLIPPING_PROFILE = "clipping_profile.py"
GENEBODY_COVERAGE = 'geneBody_coverage.py'
JUNCTION_ANNOTATION = 'junction_annotation.py'
READ_GC = 'read_GC.py'
READ_NVC = 'read_NVC.py'
READ_DISTRIBUTION = 'read_distribution.py'
READ_DUPLICATION = 'read_duplication.py'
READ_QUALITY = 'read_quality.py'

config_default = {
    'rseqc' : {
        'refgene' : config['ngs.settings']["annotation"]["transcript_annot_gtf"].replace(".gtf", ".bed12"),
        'runtime' : "01:00:00",
    },
}

update_config(config_default, config)
config = config_default


_rseqc_config_rule_default = {
    'options' : '',
    'runtime' : config['rseqc']['runtime'],
    'threads' : 1,
}
