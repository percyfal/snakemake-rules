# -*- snakemake -*-
include: "sga.settings"

config_default = {'sga' :{'preprocess' : _sga_config_rule_default.copy()}}

update_config(config_default, config)
config = config_default

rule sga_preprocess:
    """sga preprocess: filter and quality-trim reads"""
    params: cmd = config['sga']['cmd'],
            options = config['sga']['preprocess']['options'],
            runtime = config['sga']['preprocess']['runtime']
    input: reads = "{prefix}" + config['ngs.settings']['read1_label'] + config['ngs.settings']['fastq_suffix']
    output: fastq = "{prefix}.preprocess.fastq.gz",
            log = "{prefix}.preprocess.log"
    threads: config['sga']['preprocess']['threads']
    conda: "env.yaml"
    shell: "{params.cmd} preprocess {params.options} {input.reads}  2> {output.log} | gzip - > {output.fastq}"

