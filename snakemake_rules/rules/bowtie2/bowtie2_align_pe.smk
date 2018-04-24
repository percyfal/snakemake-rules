# -*- snakemake -*-
include: "bowtie2.settings.smk"

config_default = {'bowtie2' :{'align' : _bowtie2_config_rule_default.copy()}}

update_config(config_default, config)
config = config_default

rule bowtie2_align_pe:
    """Bowtie paired end alignment"""
    params: cmd = config['bowtie2']['cmd'],
            options = config['bowtie2']['align']['options'],
            index = str(config['bowtie2']['index']),
            samtools = config['samtools']['cmd'],
            runtime = config['bowtie2']['align']['runtime'],
    input: read1 = "{prefix}" + config['ngs.settings']['read1_label'] + config['ngs.settings']['fastq_suffix'],\
           read2 = "{prefix}" + config['ngs.settings']['read2_label'] + config['ngs.settings']['fastq_suffix'],
           index = expand("{index}{ext}", index=config['bowtie2']['index'], ext=config['bowtie2']['build_ext'])
    output: bam = "{prefix}.bam"
    benchmark: "{prefix}.json"
    threads: config['bowtie2']['align']['threads']
    log: "{prefix}.bwt2.log"
    conda: "env.yaml"
    shell:
        "{params.cmd} -p {threads} {params.options} -x {params.index} -1 {input.read1} -2 {input.read2} 2> {log} | {params.samtools} view -bS - > {output.bam}"
