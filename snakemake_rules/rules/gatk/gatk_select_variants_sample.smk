# -*- snakemake -*-
include: 'gatk.settings.smk'
include: "gatk_select_variants.smk"

config_default = {'gatk': {'select_variants_sample' : _gatk_config_rule_default.copy()}}

update_config(config_default, config)
config = config_default


rule gatk_select_variants_sample:
    """Run GATK SelectVariants to select variants based on a sample"""
    wildcard_constraints:
        suffix = "(.vcf|.vcf.gz)"
    params: cmd = config['gatk']['cmd'] + " -T " + SELECT_VARIANTS,
            options = " ".join(["-R", config['gatk']['select_variants']['ref'],
                                config['gatk']['select_variants']['options']]),
            runtime = config['gatk']['select_variants_sample']['runtime']
    input: vcf="{prefix}{suffix}"
    output: "{prefix}.sample_{sample}{suffix}"
    threads: config['gatk']['select_variants_sample']['threads']
    conda: "env.yaml"
    shell: "{params.cmd} {params.options} -sn {wildcards.sample} --variant {input.vcf} --out {output}"

