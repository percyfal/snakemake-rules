# -*- snakemake -*-
include: 'gatk.settings.smk'

config_default = {'gatk':{'select_variants' :_gatk_config_rule_default.copy()}}
config_default['gatk']['select_variants'].update(
    {
        'options' : " ".join(["--selectTypeToInclude SNP",
                              "--selectTypeToInclude", "INDEL",
                              "--selectTypeToInclude", "MIXED",
                              "--selectTypeToInclude", "MNP",
                              "--selectTypeToInclude", "SYMBOLIC",
                              "--selectTypeToInclude", "NO_VARIATION"]),
    })

update_config(config_default, config)
config = config_default


rule gatk_select_variants:
    """Run GATK SelectVariants to select variants"""
    wildcard_constraints:
        suffix = "(.vcf|.vcf.gz)"
    params: cmd = config['gatk']['cmd'] + " -T " + SELECT_VARIANTS,
            options = " ".join(["-R", config['gatk']['select_variants']['ref'],
                                config['gatk']['select_variants']['options']]),
            runtime = config['gatk']['select_variants']['runtime']
    input: "{prefix}{suffix}"
    output: "{prefix}.all{suffix}"
    threads: config['gatk']['select_variants']['threads']
    conda: "env.yaml"
    shell: "{params.cmd} {params.options} --variant {input} --out {output}"

