# -*- snakemake -*-
include: 'gatk.settings.smk'
include: 'gatk_variant_snp_JEXL_filtration.rule'
include: 'gatk_variant_indel_JEXL_filtration.rule'

config_default = {'gatk': {'combine_variants': _gatk_config_rule_default.copy()}}

update_config(config_default, config)
config = config_default

cmd = re.sub("-Xmx[0-9a-zA-Z]+", "-Xmx{mem}".format(mem=config['gatk']['combine_variants']['java_mem']), config['gatk']['cmd'])

rule gatk_combine_variants:
    """Run GATK CombineVariants to combine variant files.
    
    The default rule combines files with suffixes filteredSNP.vcf and
    filteredINDEL.vcf.

    """
    wildcard_constraints:
        suffix = "(.vcf|.vcf.gz)"
    params: cmd = cmd + " -T " + COMBINE_VARIANTS,
            options = " ".join(["-R", config['gatk']['combine_variants']['ref'],
                                config['gatk']['combine_variants']['options']]),
            runtime = config['gatk']['combine_variants']['runtime']
    input: "{prefix}.snp.filteredSNP{suffix}", "{prefix}.indel.filteredINDEL{suffix}"
    output: "{prefix}.variants{suffix}"
    threads: config['gatk']['combine_variants']['threads']
    conda: "env.yaml"
    shell: "command=\"{params.cmd} {params.options} $(echo {input} | sed -e 's/[^ ][^ ]*/-V &/g') -o {output}\"; eval \"${{command}}\""
