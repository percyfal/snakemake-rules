# -*- snakemake -*-
include: 'gatk.settings'

config_default = {'gatk' : {'variants_to_table' : _gatk_config_rule_default.copy()}}
config_default['gatk']['variants_to_table'].update({
    'options': "--allowMissingData",
    'fields' : ["CHROM", "POS", "FILTER", "QUAL", "DP",
                "AC", "QD", "FS", "SOR", "MQ", "MQRankSum", "ReadPosRankSum",
                "InbreedingCoeff"],
})


update_config(config_default, config)
config = config_default


rule gatk_variants_to_table:
    """Run GATK VariantToTable.

    Note: there seems to be some issue with logging. If a sequence
    dictionary is missing, a warning message is printed on stdout,
    which gets redirected to the output file. We avoid the message by
    inverted grep."""
    wildcard_constraints:
        suffix = "(.vcf|.vcf.gz)"
    params: cmd = config['gatk']['cmd'] + " -T " + VARIANTS_TO_TABLE,
            options = " ".join([
                "-R", config['gatk']['variants_to_table']['ref'],
                config['gatk']['variants_to_table']['options'], \
                " ".join("-F {}".format(x) for x in config['gatk']['variants_to_table']['fields'])]),
            runtime = config['gatk']['variants_to_table']['runtime']
    threads: config['gatk']['variants_to_table']['threads']
    conda: "env.yaml"
    input: vcf = "{prefix}{suffix}",
           vcftbi = "{prefix}{suffix}.tbi",
           ref = config['gatk']['variants_to_table']['ref'],
           refidx = config['gatk']['variants_to_table']['ref'] + ".fai",
    output: tab = "{prefix}{suffix}.table.gz"
    log: "{prefix}{suffix}.table.gz.log"
    shell: "{params.cmd} {params.options} -V {input.vcf} 2> {log} | grep -v WARN | gzip > {output.tab}"
