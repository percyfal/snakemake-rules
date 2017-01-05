# -*- snakemake -*-
include: "kpal.settings"

config_default = {'kpal' :{'distance' : _kpal_config_rule_default.copy()}}


update_config(config_default, config)
config = config_default

rule kpal_distance:
    """kpal: calculate distance between two profiles."""
    params: cmd = config['kpal']['cmd'],
            options = config['kpal']['distance']['options'],
            runtime = config['kpal']['distance']['runtime']
    wildcard_constraints: kmer_left = "[0-9]+", kmer_right = "[0-9]+"

    input: left = "{prefix_left}.k{kmer_left}", right = "{prefix_right}.k{kmer_right}"
    output: res = "{prefix_left}_{prefix_right}.txt"
    threads: config['kpal']['distance']['threads']
    conda: "env.yaml"
    shell:
        "{params.cmd} {params.options} distance {input.left} {input.right} > {output.res}"
