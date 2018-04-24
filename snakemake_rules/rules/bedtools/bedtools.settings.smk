# -*- snakemake -*-
include: '../ngs.settings.smk'

# Programs
BEDTOOLS = "bedtools"
BEDTOOLS_BAMTOBED = "bamtobed"
BEDTOOLS_COVERAGE = "coverage"
BEDTOOLS_GENOMECOV = "genomecov"
BEDTOOLS_INTERSECT = "intersect"
BEDTOOLS_MERGE = "merge"
BEDTOOLS_MAKEWINDOWS = "makewindows"

config_default = { 
    'bedtools' : {
        'cmd' : "bedtools",
        'afile' : "",
        'bfile' : "",
        'chromsizes' : "chrom.sizes",
        'options' : "",
        'sequence_capture' : {
            'bait_regions' : config['ngs.settings']['sequence_capture']['bait_regions'],
            'target_regions' : config['ngs.settings']['sequence_capture']['target_regions'],
        },
        'runtime' : "01:00:00",
    },
}

update_config(config_default, config)
config = config_default

_bedtools_config_rule_default = {
    'afile' : config['bedtools']['afile'],
    'bfile' : config['bedtools']['bfile'],
    'options' : '',
    'runtime' : config['bedtools']['runtime'],
    'threads' : 1,
}
