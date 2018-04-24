# -*- snakemake -*-
include: "bamtools.settings"

config_default = {'bamtools' :{'filter_script' : _bamtools_config_rule_default.copy()}}

update_config(config_default, config)
config = config_default


rule bamtools_filter_script:
    """Run bamtools filter on a bam file using a script input file"""
    params: cmd = config['bamtools']['cmd'],
            options = " ".join("-{} \"{}\"".format(k,v) for k,v in config['bamtools']['filter_script']['options'].items()),
            script = "-script " if config['bamtools']['filter_script']['regions'] else "",
            runtime = config['bamtools']['filter_script']['runtime']
    input: script = "{prefix}.script" if config['bamtools']['filter_script']['regions'] else [], bam = "{prefix}.bam"
    output: bam = "{prefix}.filter.bam"
    log: "{prefix}.filter.log"
    threads: config['bamtools']['filter_script']['threads']
    conda: "env.yaml"
    shell:
        "{params.cmd} filter -in {input.bam} -out {output.bam} {params.options} {params.script} {input.script} > {log}"
