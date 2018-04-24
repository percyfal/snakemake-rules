# -*- snakemake -*-
include: "bedtools.settings.smk"

config_default = {'bedtools' :{'min_genomecov' : _bedtools_config_rule_default.copy()}}

update_config(config_default, config)
config = config_default


rule bedtools_min_genomecov:
    """bedtools: extract regions with minimum coverage from input file

    The minimum coverage is defined by the mincov wildcard in the
    output file name.
    """
    params: cmd = " ".join([BEDTOOLS, BEDTOOLS_GENOMECOV]),
            mergecmd = " ".join([BEDTOOLS, BEDTOOLS_MERGE]),
            options=config['bedtools']['min_genomecov']['options'],
            runtime=config['bedtools']['min_genomecov']['runtime'],
    wildcard_constraints: suffix="(bed|vcf|gff|bam)", mincov="\d+"
    input: i = "{prefix}.{suffix}", g = config['bedtools']['chromsizes']
    output: bed = "{prefix}.{suffix}.{mincov}.coverage.bed"
    threads: config['bedtools']['min_genomecov']['threads']
    conda: "env.yaml"
    shell:
        "ext=\"\"; if [ \"{wildcards.suffix}\" = \"bam\" ]; then ext=\"bam\"; else ext=\"\"; fi\n"
        "{params.cmd} {params.options} -bga -i${{ext}} {input.i} -g {input.g} | awk -v t={wildcards.mincov} '$4>t' | {params.mergecmd} -i -  > {output.bed}"
