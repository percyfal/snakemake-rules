# -*- snakemake -*-
include: "bamtools.settings"

config_default = {'bamtools' :{'filter' : _bamtools_config_rule_default.copy()}}

update_config(config_default, config)
config = config_default


rule bamtools_filter:
    """Run bamtools filter on a bam file"""
    params: cmd = config['bamtools']['cmd'],
            options = " ".join("-{} \"{}\"".format(k,v) for k,v in config['bamtools']['filter']['options'].items()),
            runtime = config['bamtools']['filter']['runtime']
    input: bam = "{prefix}.bam"
    output: bam = "{prefix}.filter.bam"
    log: "{prefix}.filter.log"
    threads: config['bamtools']['filter']['threads']
    conda: "env.yaml"
    shell: "{params.cmd} filter -in {input.bam} -out {output.bam} {params.options} > {log}"

