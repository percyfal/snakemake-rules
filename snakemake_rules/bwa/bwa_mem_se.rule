# -*- snakemake -*-
include: "bwa.settings"
include: "bwa_index.rule"

config_default = {'bwa' :{'mem_se' : _bwa_config_rule_default.copy()}}

update_config(config_default, config)
config = config_default


rule bwa_mem_se:
    """Run bwa mem, single end"""
    params: options = config['bwa']['mem_se']['options'],
            index = config['bwa']['index'],
            runtime = config['bwa']['mem_se']['runtime']
    input: read = "{prefix}" + config['ngs.settings']['read1_label'] + config['ngs.settings']['fastq_suffix'],
           index = expand("{index}{ext}", index=config['bwa']['index'], ext=config['bwa']['index_ext'])
    output: bam = "{prefix}.bam"
    log: log = "{prefix}.bwa.log"
    threads: config['bwa']['threads']
    conda: "env.yaml"
    shell:
        "bwa mem -t {threads} {params.options} {params.index} " + \
        "{input.read} 2> {log} | " + \
        " samtools view -Sb - > {output.bam}"
