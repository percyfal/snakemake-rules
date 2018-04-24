# -*- snakemake -*-
include: "bowtie.settings.smk"

config_default = {'bowtie' :{'align' : _bowtie_config_rule_default.copy()}}
config_default['bowtie']['align'].update({'options': '--chunkmbs 200'})

update_config(config_default, config)
config = config_default

rule bowtie_align_pe:
    """Bowtie paired end alignment"""
    params: cmd = config['bowtie']['cmd'],
            options = config['bowtie']['align']['options'],
            index = config['bowtie']['index'],
            samtools = config['samtools']['cmd'],
            runtime = config['bowtie']['align']['runtime'],
    input: read1 = "{prefix}" + config['ngs.settings']['read1_label'] + config['ngs.settings']['fastq_suffix'],\
           read2 = "{prefix}" + config['ngs.settings']['read2_label'] + config['ngs.settings']['fastq_suffix'],\
           index = expand("{index}{ext}", index=config['bowtie']['index'], ext=config['bowtie']['build_ext'])
    output: bam = "{prefix}.bam"
    benchmark: "{prefix}.bwt.json"
    threads: config['bowtie']['align']['threads']
    log: "{prefix}.bwt.log"
    conda: "env.yaml"
    shell:
        #"if file --mime-type -b {input.read1} | grep -q zip; then "
        "{params.cmd} -S -p {threads} {params.options} {params.index} -1 <(gunzip -c {input.read1}) -2 <(gunzip -c {input.read2}) 2> {log} | {params.samtools} view -bS - > {output.bam}; "
        #" else "
        #t"{params.cmd} -S -p {threads} {params.options} {params.index} -1 {input.read1} -2 {input.read2} 2> {log} | {params.samtools} view -bS - > {output.bam};"
        #" fi"
