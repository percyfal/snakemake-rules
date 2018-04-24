# -*- snakemake -*-
include: "bwa.settings.smk"
include: "bwa_index.rule"
include: "bwa_aln.rule"

config_default = {'bwa' :{'samse' : _bwa_config_rule_default.copy()}}

update_config(config_default, config)
config = config_default


rule bwa_samse:
    """Run bwa samse"""
    params: options = config['bwa']['samse']['options'],
            index = config['bwa']['index'],
            runtime = config['bwa']['samse']['runtime']
    input: sai = "{prefix}.sai",
           index = config['bwa']['index'],
           read = "{prefix}" + config['ngs.settings']['read1_label'] + config['ngs.settings']['fastq_suffix']
    output: bam = "{prefix}.bam"
    log: log = "{prefix}.samse.log"
    threads: config['bwa']['samse']['threads']
    conda: "env.yaml"
    shell:
        "bwa samse {params.options} {params.index} " + \
        "{input.sai} {input.read} 2> {log} | " + \
        " samtools view -b - > {output.bam}"
