# -*- snakemake -*-
include: 'gatk.settings'
include: 'gatk_select_indel_variants.rule'

config_default = {'gatk' :{'variant_indel_JEXL_filtration' : _gatk_config_rule_default.copy()}}
config_default['gatk']['variant_indel_JEXL_filtration'].update(
    {
        'expressions':
        ["QD < 2.0", "ReadPosRankSum < -20.0", "FS > 200.0",
         "SOR > 10.0"]
    })

update_config(config_default, config)
config = config_default


rule gatk_variant_indel_JEXL_filtration:
    """Run GATK VariantFiltration

    Perform hard filtering using JEXL expressions
    """
    wildcard_constraints:
        suffix = "(.vcf|.vcf.gz)"
    params: cmd = config['gatk']['cmd'] + " -T " + VARIANT_FILTRATION,
            options = " ".join([
                " ".join(["--filterName GATKStandard{e} --filterExpression \"{exp}\"".format(e=exp.split()[0], exp=exp) \
                                for exp in config['gatk']['variant_indel_JEXL_filtration']['expressions']])
                ]),
            quote = "" if config['gatk']['cmd'].startswith("gatk") else "",
            runtime = config['gatk']['variant_indel_JEXL_filtration']['runtime']
    input: vcf = "{prefix}{suffix}", ref = config['gatk']['variant_indel_JEXL_filtration']['ref']
    output: vcf = "{prefix}.filteredINDEL{suffix}"
    threads: config['gatk']['variant_indel_JEXL_filtration']['threads']
    conda: "env.yaml"
    shell: "{params.cmd} {params.quote}{params.options}{params.quote} -R {input.ref} --variant {input.vcf} --out {output.vcf}"
