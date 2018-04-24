# -*- snakemake -*-
include: 'gatk.settings.smk'

config_default = {'gatk' :{'select_indel_variants' : _gatk_config_rule_default.copy()}}
config_default['gatk']['select_indel_variants'].update(
    {
        'options' : " ".join(["--selectTypeToInclude", "INDEL",
                              "--selectTypeToInclude", "MIXED",
                              "--selectTypeToInclude", "MNP",
                              "--selectTypeToInclude", "SYMBOLIC",
                              "--selectTypeToInclude", "NO_VARIATION"]),
    })


update_config(config_default, config)
config = config_default


rule gatk_select_indel_variants:
    """Run GATK SelectVariants to select INDELs"""
    wildcard_constraints:
        suffix = "(.vcf|.vcf.gz)"
    params: cmd = config['gatk']['cmd'] + " -T " + SELECT_VARIANTS,
            options =  config['gatk']['select_indel_variants']['options'],
            runtime =  config['gatk']['select_indel_variants']['runtime']
    input: vcf = "{prefix}{suffix}", ref = config['gatk']['select_indel_variants']['ref']
    output: vcf = "{prefix}.indel{suffix}"
    threads: config['gatk']['select_indel_variants']['threads']
    conda: "env.yaml"
    shell: "{params.cmd} {params.options} -R {input.ref} --variant {input.vcf} --out {output.vcf}"

