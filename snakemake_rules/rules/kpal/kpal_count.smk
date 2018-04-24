# -*- snakemake -*-
include: "kpal.settings.smk"

config_default = {'kpal' :{'count' : _kpal_config_rule_default.copy()}}

update_config(config_default, config)
config = config_default

rule kpal_count:
    """kpal: count k-mers."""
    params: cmd = config['kpal']['cmd'],
            options = config['kpal']['count']['options'],
            runtime = config['kpal']['count']['runtime']
    wildcard_constraints: read = config['ngs.settings']['read_label_re'],
                          kmer = "[0-9]+"
    input: read = "{prefix}{read}" + config['ngs.settings']['fastq_suffix']
    output: kmer = "{prefix}{read}.k{kmer}"
    threads: config['kpal']['count']['threads']
    conda: "env.yaml"
    shell:
        "if file --mime-type -b -L {input.read} | grep -q gzip; then zcat {input.read} | {params.cmd} count {params.options} -k {wildcards.kmer} - {output.kmer};"
        "else "
        "{params.cmd} count {params.options} -k {wildcards.kmer} {input.read} {output.kmer}; fi;"
