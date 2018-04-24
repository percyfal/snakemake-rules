# -*- snakemake -*-
include: "sga.settings.smk"
include: "../seqtk/seqtk.settings.smk"

config_default = {'sga' :{'preprocess' : _sga_config_rule_default.copy()}}

update_config(config_default, config)
config = config_default

rule sga_preprocess_pe:
    """sga preprocess: filter and quality-trim reads, paired end reads"""
    params: cmd = config['sga']['cmd'],
            seqtk = config['seqtk']['cmd'],
            options = config['sga']['preprocess']['options'],
            runtime = config['sga']['preprocess']['runtime'],
            tmp = "sga.preprocess.tmp.fastq.gz"
    input: read1 = "{prefix}" + config['ngs.settings']['read1_label'] + config['ngs.settings']['fastq_suffix'],
           read2 = "{prefix}" + config['ngs.settings']['read2_label'] + config['ngs.settings']['fastq_suffix'],
    output: fastq1 = "{prefix}" + config['ngs.settings']['read1_label'] + ".preprocess" + config['ngs.settings']['fastq_suffix'],
            fastq2 = "{prefix}" + config['ngs.settings']['read2_label'] + ".preprocess" + config['ngs.settings']['fastq_suffix'],
            log = "{prefix}.preprocess.log"
    threads: config['sga']['preprocess']['threads']
    conda: "env.yaml"
    shell: "{params.cmd} preprocess --pe-mode 1 {params.options} {input.read1} {input.read2} 2> {output.log} | gzip - > {params.tmp}; "
           "{params.seqtk} seq -1 {params.tmp} | gzip - > {output.fastq1}; "
           "{params.seqtk} seq -2 {params.tmp} | gzip - > {output.fastq2}; "
           "rm -f {params.tmp}"
           

