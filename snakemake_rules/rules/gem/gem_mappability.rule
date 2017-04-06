# -*- snakemake -*-
include: "gem.settings"

config_default = {'gem' :{'mappability' : _gem_config_rule_default.copy()}}
config_default['gem']['mappability'].update({'options' : "-l 100 -T 7", "cmd": 'gem-mappability'})

update_config(config_default, config)
config = config_default

rule gem_mappability:
    params: cmd = config['gem']['mappability']['cmd'],
            options = config['gem']['mappability']['options'],
            runtime = config['gem']['mappability']['runtime']
    input: gem = "{prefix}.gem"
    output: mappability = "{prefix}.mappability"
    threads: config['gem']['mappability']['threads']
    shell: "{params.cmd} {params.options} -I {input.gem} -o {output.mappability}"
