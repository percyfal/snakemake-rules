# -*- snakemake -*-
include: '../ngs.settings.smk'


config_default = { 
    'cutadapt' : {
        'cmd' : "cutadapt",
        'threeprime': "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC",
        'fiveprime' : "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT",
        'runtime' : '24:00:00',
        'options' : "",
        'fastq_suffix' : config['ngs.settings']['fastq_suffix'],
        'fastq_suffix_re' : config['ngs.settings']['fastq_suffix_re'],
        'read_label_re' : config['ngs.settings']['read_label_re'],
        'read1_label_re' : config['ngs.settings']['read1_label_re'],
        'read2_label_re' : config['ngs.settings']['read2_label_re'],
        'read1_label' : config['ngs.settings']['read1_label'],
        'read2_label' : config['ngs.settings']['read2_label'],
        'threads' : 1,
    },
}

update_config(config_default, config)
config = config_default


_cutadapt_config_rule_default = {
    'options' : '',
    'runtime' : config['cutadapt']['runtime'],
    'threads' : 1,
}
