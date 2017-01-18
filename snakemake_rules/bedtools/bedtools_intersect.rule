# -*- snakemake -*-
include: "bedtools.settings"

config_default = {'bedtools' :{'intersect' : _bedtools_config_rule_default.copy()}}

update_config(config_default, config)
config = config_default


rule bedtools_intersect_make_region_baits:
    """Generate overlapping bed file for bait definition file"""
    params: cmd = " ".join([BEDTOOLS, BEDTOOLS_INTERSECT]),
            options = config['bedtools']['intersect']['options'],
            runtime = config['bedtools']['intersect']['runtime'],
    input: a="{prefix}.region_{gene}.bed", b=config['bedtools']['sequence_capture']['bait_regions']
    output: "{prefix}.region_{gene}.baits.bed"
    threads: config['bedtools']['intersect']['threads']
    conda: "env.yaml"
    shell: "{params.cmd} -a {input.a} -b {input.b} > {output}"

rule bedtools_intersect_make_region_targets:
    """Generate overlapping bed file for target definition file"""
    params: cmd = " ".join([BEDTOOLS, BEDTOOLS_INTERSECT]),
            options = config['bedtools']['intersect']['options'],
            runtime = config['bedtools']['intersect']['runtime'],
    input: a="{prefix}.region_{gene}.bed", b=config['bedtools']['sequence_capture']['target_regions']
    output: "{prefix}.region_{gene}.targets.bed"
    threads: config['bedtools']['intersect']['threads']
    conda: "env.yaml"
    shell: "{params.cmd} -a {input.a} -b {input.b} > {output}"

