# -*- snakemake -*-
include: 'gatk.settings'

config_default = {'gatk': {'select_snp_variants': _gatk_config_rule_default.copy()}}
config_default['gatk']['select_snp_variants'].update({'options' : "--selectTypeToInclude SNP"})

update_config(config_default, config)
config = config_default


rule gatk_select_snp_variants:
    """Run GATK SelectVariants to select SNPs"""
    wildcard_constraints:
        suffix = "(.vcf|.vcf.gz)"
    params: cmd = config['gatk']['cmd'] + " -T " + SELECT_VARIANTS,
            options = config['gatk']['select_snp_variants']['options'],
            runtime = config['gatk']['select_snp_variants']['runtime']
    input: vcf = "{prefix}{suffix}", ref = config['gatk']['select_snp_variants']['ref']
    output: vcf = "{prefix}.snp{suffix}"
    threads: config['gatk']['select_snp_variants']['threads']
    conda: "env.yaml"
    shell: "{params.cmd} {params.options} -R {input.ref} --variant {input.vcf} --out {output.vcf}"

