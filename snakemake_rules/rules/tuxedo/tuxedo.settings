# -*- snakemake -*-
include: '../ngs.settings'

config_default = { 
    "tuxedo" : {
        'ref' : config['ngs.settings']['db']['ref'],
        'rg_fn' : None,
        'version2' : True, # Run version 2 by default
        'index' : "",
        'cufflinks' : {
            'cmd' : 'cufflinks', 
            'options' : '',
            'transcript_annot_gtf' : config['ngs.settings']['annotation']['transcript_annot_gtf'],
            'annot_label' : config['ngs.settings']['db']['build'],
            'output_dir' : os.curdir,
        },
        'build' : {
            'ext_v1' : [".1.ebwt", ".2.ebwt", ".3.ebwt", ".4.ebwt", ".rev.1.ebwt", ".rev.2.ebwt"],
            'ext_v2' : [".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2"],
        },
    },
}

update_config(config_default, config)
config = config_default

bowtie = 'bowtie2' if config['tuxedo']['version2'] else 'bowtie'
