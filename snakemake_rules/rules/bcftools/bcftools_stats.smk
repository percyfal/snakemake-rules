# -*- snakemake -*-
include: "bcftools.settings.smk"
include: "../htslib/htslib_bgzip.rule"

config_default = {'bcftools' :{'stats' : _bcftools_config_rule_default.copy()}}

update_config(config_default, config)
config = config_default

rule bcftools_stats:
    """bcftools: calculate statistics.

    NB: previously called 'vcfcheck'"""
    params: cmd = config['bcftools']['cmd'],
            options = config['bcftools']['stats']['options'],
            runtime = config['bcftools']['stats']['runtime']
    wildcard_constraints: suffix = "(vcf|vcf.gz|fofn)"
    input: vcf = "{prefix}.{suffix}",
           ref = config['bcftools']['ref'] if  config['bcftools']['ref'] else []
    output: stats = "{prefix}.{suffix}.stats"
    threads: config['bcftools']['stats']['threads']
    conda: "env.yaml"
    shell:
        "ref=\"\"; if [[ {input.ref} != \"\" ]]; then ref=\"-F {input.ref}\"; fi;\n"
        "if [[ {wildcards.suffix} == \"fofn\" ]]; then\n"
        "command=\"{params.cmd} stats {params.options} ${{ref}} $(cat {input.vcf} | tr \"\\n\", \" \") > {output.stats}\"; eval ${{command}};\n"
        "else\n"
        "{params.cmd} stats {params.options} ${{ref}} {input.vcf} > {output.stats}\n"
        "fi"

