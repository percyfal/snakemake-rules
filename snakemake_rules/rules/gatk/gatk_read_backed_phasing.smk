# -*- snakemake -*-
include: 'gatk.settings.smk'

def _gatk_read_backed_phasing_inputs_fn(wildcards):
    raise RuleException(
        "Implement function for generating source files and set the configuration parameter config['gatk']['read_backed_phasing']['inputfun'] to the name of the function. The function should return a vcf input file, a bam file and the corresponding bam index")

config_default = {'gatk': {'read_backed_phasing' : _gatk_config_rule_default.copy()}}
config_default['gatk']['read_backed_phasing'].update({'inputfun' : _gatk_read_backed_phasing_inputs_fn})

update_config(config_default, config)
config = config_default


rule gatk_read_backed_phasing:
    """Run GATK ReadBackedPhasing"""
    wildcard_constraints:
        suffix = "(.vcf|.vcf.gz)"
    params: cmd = config['gatk']['cmd'] + " -T " + READ_BACKED_PHASING,
            options = " ".join(["-R", config['gatk']['read_backed_phasing']['ref'],
                                config['gatk']['read_backed_phasing']['options']]),
            runtime = config['gatk']['read_backed_phasing']['runtime'],
    input: config['gatk']['read_backed_phasing']['inputfun']
    output: "{prefix}.phased{suffix}"
    threads: config['gatk']['read_backed_phasing']['threads']
    conda: "env.yaml"
    shell: "{params.cmd} {params.options} -I {input[1]} -o {output} --variant {input[0]}"
