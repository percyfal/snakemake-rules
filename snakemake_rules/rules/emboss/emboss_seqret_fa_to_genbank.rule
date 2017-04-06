# -*- snakemake -*-
include: "emboss.settings"

config_default = {'emboss' :{'seqret_fa_to_genbank' : _emboss_config_rule_default.copy()}}

update_config(config_default, config)
config = config_default


rule emboss_seqret_fa_to_genbank:
    """EMBOSS seqret: convert fasta, suffix .fa, to genbank"""
    params: runtime = config['emboss']['seqret_fa_to_genbank']['runtime']
    input: "{prefix}.fa"
    output: "{prefix}.genbank"
    conda: "env.yaml"
    threads: config['emboss']['seqret_fa_to_genbank']['threads']
    shell: "seqret -sequence fasta::{input} -outseq genbank::{output}"

