# -*- snakemake -*-
include: '../ngs.settings'

config_default = { 
    'rsem' : {
        'annot_gtf' : config['ngs.settings']['annotation']['transcript_annot_gtf'],
        'threads' : config['settings']['threads'],
        'ref' : config['ngs.settings']['db']['ref'],
        'ref_sfx' : '.transcripts.fa',
        'index' : "rsem_index",
        'index_is_transcriptome': False, # Set to true if mapping directly to transcriptome
        'use_multimapped': True, # RSEM uses multimapped reads for quantification
        'calculate-expression' : {
            'cmd' : 'rsem-calculate-expression',
            'options' : '--no-bam-output',
            'bowtie-options' : "--no-bam-output --bowtie-chunkmbs 512",
        },
        'runtime' : "01:00:00",
    },
}

update_config(config_default, config)
config = config_default

# Add rsem to rnaseq section
config['ngs.settings']['rnaseq']['_quantification'].append('rsem')
# Set transcriptome based quantification true (for STAR etc)
if not config['rsem']['index_is_transcriptome']:
    config['ngs.settings']['rnaseq']['_transcriptome_quantification'] = True

