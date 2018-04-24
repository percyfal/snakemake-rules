# -*- snakemake -*-
include: "vcflib.settings.smk"

config_default = {'vcflib' :{'check' : _vcflib_config_rule_default.copy()}}

update_config(config_default, config)
config = config_default


rule vcfcheck:
    """vcfcheck: Run vcfcheck on vcf file"""
    params: cmd='vcfcheck',
            runtime=config['vcflib']['check']['runtime']
    input: vcf = "{prefix}.vcf",
           ref = config['vcflib']['ref']
    output: "{prefix}.check"
    threads: config['vcflib']['check']['threads']
    conda: "env.yaml"
    shell: "{params.cmd} {input.vcf} -f {input.ref} > {output}"
