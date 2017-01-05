# -*- snakemake -*-
include: '../ngs.settings'

config_default = {
    'rpkmforgenes' : {
        'cmd' : 'rpkmforgenes.py',
        'annotation' : config['ngs.settings']['annotation']['transcript_annot_gtf'],
        'options' : "-readcount -fulltranscript -mRNAnorm -rmnameoverlap -bothendsceil",
        'unique' : None,
        'use_multimapped': False, # rpkmforgenes requires uniquely mapping reads
        'annot_format' : 'refFlat',
        'runtime' : '01:00:00',
        'threads' : 1,
    },
}

update_config(config_default, config)
config = config_default

# Add rpkmforgenes to rnaseq section
config['ngs.settings']['rnaseq']['_quantification'].append('rpkmforgenes')
