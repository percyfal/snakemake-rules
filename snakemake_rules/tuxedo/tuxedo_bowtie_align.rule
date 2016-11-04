# -*- snakemake -*-
include: "tuxedo.settings"
include: "../samtools/samtools.settings"

config_default = {
    'tuxedo' : {
        'bowtie' : {
            'threads' : config['settings']['threads'],
            'options' : '',
        },
    },
}

update_config(config_default, config)
config = config_default


rule tuxedo_bowtie_align:
    """Bowtie paired end alignment"""
    params: cmd = bowtie,
            options = config['tuxedo']['bowtie']['options'],
            index = ("-x " if config['tuxedo']['version2'] else "") + str(config['tuxedo']['index'])
    input: read1="{prefix}" + config['ngs.settings']['read1_label'] + config['ngs.settings']['fastq_suffix'],\
           read2="{prefix}" + config['ngs.settings']['read2_label'] + config['ngs.settings']['fastq_suffix']
    output: "{prefix}.bam"
    threads: config['tuxedo']['bowtie']['threads']
    shell: "{params.cmd} -p {threads} {params.options} {params.index} -1 {input.read1} -2 {input.read2} | " + config['samtools']['cmd'] + " view -Sb - > {output}.tmp && mv {output}.tmp {output}"

