# -*- snakemake -*-
include: 'gatk.settings.smk'
include: '../bcftools/bcftools_stats.smk'

config_default = {'gatk' : {'vqsr_prepare_training_data' :  _gatk_config_rule_default.copy()}}
config_default['gatk']['vqsr_prepare_training_data'].update({'fraction': 0.2})

update_config(config_default, config)
config = config_default


rule gatk_vqsr_prepare_training_data:
     """Prepare VQSR training data by selecting variants in top fraction based on quality.

    Fraction is set in
    config['gatk']['vqsr_prepare_training_data']['fraction'] and
    defaults to 0.2.

    """
    wildcard_constraints:
        mode = "(snp|indel)"
    params: options = config['gatk']['vqsr_prepare_training_data']['options'],
            runtime = config['gatk']['vqsr_prepare_training_data']['runtime'],
            fraction = config['gatk']['vqsr_prepare_training_data']['fraction'],
            bcftools = config['bcftools']['cmd']
    input: vcf = "{prefix}.{mode}.vcf.gz", stats = "{prefix}.{mode}.vcf.gz.stats"
    output: "{prefix}.{mode}.train.vcf.gz"
    threads: config['gatk']['vqsr_prepare_training_data']['threads']
    log: "{prefix}.{mode}.train.log"
    shell:
        "N=$(grep records {input.stats} | cut -f 4);\n"
        "pos=$(echo \"$N * {params.fraction} / 1\" | bc);\n"
        "qual=$({params.bcftools} query -f '%QUAL\\n' {input.vcf} | sort -g -r | head -$pos | tail -1);\n"
        "{params.bcftools} view -i \"MIN(QUAL>$qual)\" {input.vcf} -O z -o {output}"
