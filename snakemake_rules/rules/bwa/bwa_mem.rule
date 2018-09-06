# -*- snakemake -*-
include: "bwa.settings"
include: "bwa_index.rule"

config_default = {'bwa' :{'mem' : _bwa_config_rule_default.copy()}}

update_config(config_default, config)
config = config_default

rule bwa_mem:
    """Run bwa mem"""
    params: options = config['bwa']['mem']['options'],
            index = config['bwa']['index'],
            runtime = config['bwa']['mem']['runtime']
    input: read1 = "{prefix}" + config['ngs.settings']['read1_label'] + config['ngs.settings']['fastq_suffix'],
           read2 = "{prefix}" + config['ngs.settings']['read2_label'] + config['ngs.settings']['fastq_suffix'],
           index = expand("{index}{ext}", index=config['bwa']['index'], ext=config['bwa']['index_ext'])
    output: bam = "{prefix}.bam"
    log: log = "{prefix}.log"
    threads: config['bwa']['mem']['threads']
    conda: "env.yaml"
    shell:
        "bwa mem -t {threads} {params.options} {params.index} " + \
        "{input.read1} {input.read2} 2> {log} | " + \
        " samtools view -Sb - > {output.bam}"
