# -*- snakemake -*-
include: "seqtk.settings"

config_default = {'seqtk' :{'seq' : _seqtk_config_rule_default.copy()}}

update_config(config_default, config)
config = config_default

rule seqtk_seq_fastq_to_fasta:
    """seqtk: convert fastq to fasta"""
    params: cmd = config['seqtk']['cmd'],
            options = config['seqtk']['seq']['options'],
            runtime = config['seqtk']['seq']['runtime']
    input: reads = "{prefix}" + config['ngs.settings']['fastq_suffix']
    output: fasta = "{prefix}.fasta"
    threads: config['seqtk']['seq']['threads']
    conda: "env.yaml"
    shell: "{params.cmd} seq {params.options} -a {input.reads} > {output.fasta}"

