# -*- snakemake -*-
include: "kpal.settings.smk"

config_default = {'kpal' :{'count_all' : _kpal_config_rule_default.copy()}}
config_default['kpal']['count_all'].update({'input' : []})

update_config(config_default, config)
config = config_default

rule kpal_count_all:
    """kpal: count k-mers for multiple input files"""
    params: cmd = config['kpal']['cmd'],
            options = config['kpal']['count_all']['options'],
            runtime = config['kpal']['count_all']['runtime']
    wildcard_constraints: kmer = "[0-9]+"
    input: config['kpal']['count_all']['input']
    output: kmer = os.path.join("{path}", "kpal{label}.k{kmer}")
    threads: config['kpal']['count_all']['threads']
    conda: "env.yaml"
    shell:
        "command=\"{params.cmd} count {params.options} -k {wildcards.kmer} $(echo {input} | sed -e 's/[^ ][^ ]*/<(seqtk seq -a &)/g') {output.kmer}\"; eval \"${{command}}\""
