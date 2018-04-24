# -*- snakemake -*-
include: "bwa.settings.smk"
include: "bwa_index.rule"

config_default = {'bwa' :{'aln' : _bwa_config_rule_default.copy()}}

update_config(config_default, config)
config = config_default


rule bwa_aln:
    """Run bwa aln"""
    params: options = config['bwa']['aln']['options'],
            index = config['bwa']['index'],
            runtime = config['bwa']['aln']['runtime']
    input: read = "{prefix}" + config['ngs.settings']['read1_label'] + config['ngs.settings']['fastq_suffix'],
           index = expand("{index}{ext}", index=config['bwa']['index'], ext=config['bwa']['index_ext'])
    output: bam = "{prefix}.sai"
    log: log = "{prefix}.bwa.log"
    threads: config['bwa']['aln']['threads']
    conda: "env.yaml"
    shell:
        "bwa aln -t {threads} {params.options} {params.index} " + \
        "{input.read} -f {output} 2> {log}"
