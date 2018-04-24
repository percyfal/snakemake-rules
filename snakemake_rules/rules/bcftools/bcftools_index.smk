# -*- snakemake -*-
include: "bcftools.settings.smk"

config_default = {'bcftools' :{'index' : _bcftools_config_rule_default.copy()}}

update_config(config_default, config)
config = config_default

rule bcftools_index:
    """bcftools index: index variants"""
    wildcard_constraints:
        sfx = "(.vcf.gz|.bcf)"
    params: cmd = config['bcftools']['cmd'],
            options = config['bcftools']['index']['options'],
            runtime = config['bcftools']['index']['runtime']
    input: vcf = "{prefix}{sfx}"
    output: csi = "{prefix}{sfx}.csi"
    threads: config['bcftools']['index']['threads']
    conda: "env.yaml"
    shell: "{params.cmd} index {params.options} {input.vcf}"

