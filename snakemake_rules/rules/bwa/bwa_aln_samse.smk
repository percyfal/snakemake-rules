# -*- snakemake -*-
include: "bwa.settings.smk"
include: "bwa_index.smk"

config_default = {'bwa' : {'aln' : _bwa_config_rule_default.copy(),
                           'samse' : _bwa_config_rule_default.copy()}}


update_config(config_default, config)
config = config_default


rule bwa_aln_samse:
    """Run bwa aln followed by samse; outputs bam"""
    params: options = config['bwa']['aln']['options'],
            index = config['bwa']['index'],
            runtime = config['bwa']['aln']['runtime'],
            options_samse = config['bwa']['samse']['options']
    input: read = "{prefix}" + config['ngs.settings']['read1_label'] + config['ngs.settings']['fastq_suffix'],
           index = expand("{index}{ext}", index=config['bwa']['index'], ext=config['bwa']['index_ext'])
    output: bam = "{prefix}.bam"
    log: log = "{prefix}.bwa.log"
    threads: config['bwa']['aln']['threads']
    conda: "env.yaml"
    shell:
        "bwa aln -t {threads} {params.options} {params.index} " + \
        "{input.read} 2> {log} | bwa samse {params.options_samse} {params.index} - {input.read}" + \
        "| samtools view -b - > {output.bam}"
