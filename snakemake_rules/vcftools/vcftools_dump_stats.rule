# -*- snakemake -*-
include: "vcftools.settings"

config_default = {'vcftools' :{'dump_stats' : _vcftools_config_rule_default.copy()}}

update_config(config_default, config)
config = config_default

rule vcftools_dump_stats:
    """vcftools: dump statistics"""
    params: cmd='vcf-stats',
            options=config['vcftools']['dump_stats']['options'],
            runtime = config['vcftools']['dump_stats']['runtime']
    input: "{prefix}.vcf"
    output: "{prefix}.stats.dump"
    threads: config['vcftools']['dump_stats']['threads']
    shell: "{params.cmd} {params.options} {input} -p {wildcards.prefix}.stats"

