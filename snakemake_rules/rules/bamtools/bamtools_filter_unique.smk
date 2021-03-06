# -*- snakemake -*-
include: "bamtools.settings.smk"

config_default = {'bamtools' :{'filter_unique' : _bamtools_config_rule_default.copy()}}

update_config(config_default, config)
config = config_default


rule bamtools_filter_unique:
    """Run bamtools filter on a bam file"""
    params: cmd = config['bamtools']['cmd'],
            options = " ".join("-{} \"{}\"".format(k,v) for k,v in config['bamtools']['filter_unique']['options'].items()),
            runtime = config['bamtools']['filter_unique']['runtime']
    input: bam = "{prefix}.bam"
    output: bam = "{prefix}_unique.bam"
    log: "{prefix}_unique.log"
    threads: config['bamtools']['filter_unique']['threads']
    conda: "env.yaml"
    shell: "{params.cmd} filter -in {input.bam} -out {output.bam} {params.options} > {log}"

