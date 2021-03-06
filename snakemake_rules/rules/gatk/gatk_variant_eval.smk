# -*- snakemake -*-
include: 'gatk.settings.smk'

config_default = {'gatk' :{'variant_eval': _gatk_config_rule_default.copy()}}
config_default['gatk']['variant_eval'].update(
    {
        'options' : " ".join(["-ST Filter -l INFO --doNotUseAllStandardModules --evalModule CompOverlap --evalModule CountVariants --evalModule TiTvVariantEvaluator --evalModule ValidationReport"])
    })

update_config(config_default, config)
config = config_default


rule gatk_variant_eval:
    """Run GATK VariantEval"""
    wildcard_constraints:
        suffix = "(.vcf|.vcf.gz)"
    params: cmd = config['gatk']['cmd'] + " -T " + VARIANT_EVAL,
            options = " ".join(["-R", config['gatk']['variant_eval']['ref'],
                                config['gatk']['variant_eval']['options'],
                                "--dbsnp {known}".format(known=config['gatk']['known_sites'] if not config['gatk']['known_sites'] == "" else "")]),
            runtime = config['gatk']['variant_eval']['runtime']
    input: "{prefix}{suffix}"
    output: "{prefix}{suffix}.eval_metrics"
    threads: config['gatk']['variant_eval']['threads']
    conda: "env.yaml"
    shell: "{params.cmd} {params.options} --eval {input} -o {output}"

