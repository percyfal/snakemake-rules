# -*- snakemake -*-
include: 'gatk.settings.smk'
include: 'gatk_select_snp_variants.smk'

config_default = {'gatk' :{'variant_snp_JEXL_filtration' : _gatk_config_rule_default.copy()}}
config_default['gatk']['variant_snp_JEXL_filtration'].update(
    {
        'expressions': ["QD < 2.0", "MQ < 40.0", "FS > 60.0",
                        "MQRankSum < -12.5",
                        "ReadPosRankSum < -8.0", "SOR > 3.0"],
    })


update_config(config_default, config)
config = config_default


rule gatk_variant_snp_JEXL_filtration:
    """Run GATK VariantFiltration on SNPs

    Perform hard filtering using JEXL expressions.
    """
    wildcard_constraints:
        suffix = "(.vcf|.vcf.gz)"
    params: cmd = config['gatk']['cmd'] + " -T " + VARIANT_FILTRATION,
            options = " ".join([
                " ".join(["--filterName GATKStandard{e} --filterExpression \"{exp}\"".format(e=exp.split()[0], exp=exp) \
                                for exp in config['gatk']['variant_snp_JEXL_filtration']['expressions']])
                ]),
            quote = "" if config['gatk']['cmd'].startswith("gatk") else "",
            runtime = config['gatk']['variant_snp_JEXL_filtration']['runtime']
    input: vcf = "{prefix}{suffix}", ref = config['gatk']['variant_snp_JEXL_filtration']['ref']
    output: vcf = "{prefix}.filteredSNP{suffix}"
    threads: config['gatk']['variant_snp_JEXL_filtration']['threads']
    conda: "env.yaml"
    shell: "{params.cmd} {params.quote}{params.options}{params.quote} -R {input.ref} --variant {input.vcf} --out {output.vcf}"

