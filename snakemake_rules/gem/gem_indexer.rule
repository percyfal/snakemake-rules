# -*- snakemake -*-
include: "gem.settings"

config_default = {'gem' :{'indexer' : _gem_config_rule_default.copy()}}
config_default['gem']['indexer'].update({'cmd': 'gem-indexer', 'options' : "-T 7"})

update_config(config_default, config)
config = config_default

rule gem_indexer:
    params: cmd = config['gem']['indexer']['cmd'],
            options = config['gem']['indexer']['options'],
            runtime = config['gem']['indexer']['runtime'],
    input: fasta = "{prefix}.fasta"
    output: gem = "{prefix}.gem"
    log: log = "{prefix}.gem.log"
    threads: config['gem']['indexer']['threads']
    shell: "{params.cmd} {params.options} -i {input.fasta} -o {output.gem} > {log.log}"
