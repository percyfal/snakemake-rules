# -*- snakemake -*-
include: "emboss.settings.smk"

config_default = {'emboss' :{'seqret_fasta_to_genbank' : _emboss_config_rule_default.copy()}}

update_config(config_default, config)
config = config_default


rule emboss_seqret_fasta_to_genbank:
    """EMBOSS seqret: convert fasta to genbank"""
    params: runtime = config['emboss']['seqret_fasta_to_genbank']['runtime']
    input: "{prefix}.fasta"
    output: "{prefix}.genbank"
    conda: "env.yaml"
    threads: config['emboss']['seqret_fasta_to_genbank']['threads']
    shell: "seqret -sequence {input} -outseq genbank:{output}"
    
