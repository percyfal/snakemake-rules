# -*- snakemake -*-
include: 'gatk.settings.smk'

config_default = {'gatk' :{'variant_recalibrator': _gatk_config_rule_default.copy()}}
config_default['gatk']['variant_recalibrator'].update(
    {
        'options': "-allPoly -mG 4 -std 10.0",
        'resources': {'snp': {}, 'indel': {}},
        'annotations': ["DP", "QD", "FS", "SOR", "MQ", "MQRankSum", "ReadPosRankSum", "InbreedingCoeff"],
        'tranches': [100.0, 99.9, 99.0, 90.0]
    })

update_config(config_default, config)
config = config_default


rule gatk_variant_recalibrator:
    """Run GATK VariantRecalibrator.

    This rule contains a lot of options. For clarity, the options
    configuration has been partitioned into four different
    configuration settings:

    options (str): regular options
    annotations (list): list of annotations to use
    tranches (list): list of tranches

    The run command sets resources according to the following option:
    resources (dict): dictionary of resource mode:{file:option} pairs where
    the options specify truth/training/known status along with a label
    and prior

    """
    wildcard_constraints:
        mode = "(snp|indel)"
    params: cmd = config['gatk']['cmd'] + " -T " + VARIANT_RECALIBRATOR,
            options =  config['gatk']['variant_recalibrator']['options'],
            annotations = " ".join("-an {}".format(x) for x in config['gatk']['variant_recalibrator']['annotations']),
            tranches = " ".join("-tranche {}".format(x) for x in config['gatk']['variant_recalibrator']['tranches']),
            runtime = config['gatk']['variant_recalibrator']['runtime']
    input: vcf = "{prefix}.vcf.gz", tbi = "{prefix}.vcf.gz.tbi",
           ref = config['gatk']['variant_recalibrator']['ref'],
           fai = config['gatk']['variant_recalibrator']['ref'] + ".fai",
           resources = [x for k,v in config['gatk']['variant_recalibrator']['resources'].items() for x in v.keys()],
           resources_tbi = ["{}.tbi".format(x) for k,v in config['gatk']['variant_recalibrator']['resources'].items() for x in v.keys()]
    output: recal="{prefix}.{mode}.recal", tranchesFile="{prefix}.{mode}.tranches", rscript="{prefix}.{mode}.R"
    threads: config['gatk']['variant_recalibrator']['threads']
    run:
        resources = " ".join("-resource:{} {}".format(v, k) for k,v in config['gatk']['variant_recalibrator']['resources'][wildcards.mode].items() if k.startswith(wildcards.prefix))
        d = dict(params)
        d.update(dict(input))
        d.update(dict(output))
        d.update(dict(wildcards))
        d.update({'resources': resources,
                  'threads': threads})
        d['mode'] = d['mode'].upper()
        shell("{cmd} -R {ref} -nt {threads} -mode {mode} -input {vcf} {resources} {annotations} {tranches} {options} -recalFile {recal} -tranchesFile {tranchesFile} -rscriptFile {rscript}".format(**d))
