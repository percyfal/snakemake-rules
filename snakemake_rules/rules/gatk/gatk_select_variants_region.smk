# -*- snakemake -*-
include: 'gatk.settings.smk'
include: "gatk_select_variants.smk"
include: "gatk_unified_genotyper.smk"

config_default = {'gatk' :{'select_variants_region' : _gatk_config_rule_default.copy()}}

update_config(config_default, config)
config = config_default


rule gatk_select_variants_region:
    """Run GATK SelectVariants to select variants based on a region"""
    wildcard_constraints:
        suffix = "(.vcf|.vcf.gz)"
    params: cmd = config['gatk']['cmd'] + " -T " + SELECT_VARIANTS,
            options = " ".join(["-R", config['gatk']['select_variants']['ref'],
                                config['gatk']['select_variants']['options']]),
            runtime = config['gatk']['select_variants_region']['runtime']
    input: vcf="{prefix}{suffix}", bed="{prefix}.region_{region}.bed"
    output: "{prefix}.region_{region}{suffix}"
    threads: config['gatk']['select_variants_region']['threads']
    conda: "env.yaml"
    shell: "{params.cmd} {params.options} -L {input.bed} --variant {input.vcf} --out {output}"

ruleorder: gatk_select_variants_region > gatk_unified_genotyper
