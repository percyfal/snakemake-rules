# -*- snakemake -*-
#
# Manual installation required. No setup exists so full path to
# installation directory must be supplied. Download from
# https://sites.google.com/site/danposdoc/download and extract archive
#
# Requirements:
# R version 2.9.1.
# Python 2.7, rpy2, numpy 1.5.0.
# Samtools 0.1.7 (only when input data is in sam or bam format)
#
# Memory ~ (genome_size/step_size) x ( replicate_count + max(2,
# comparison_count) ) x 8 bits.
#
include: "../ngs.settings.smk"


config_default = {
    'danpos' : {
        'bins': [(0,100), (180,247), (315,473), (558,615)],
        'cmd': 'danpos.py',
    },
}

update_config(config_default, config)
config = config_default
