# -*- snakemake -*-
include: "htslib.settings"

config_default = {'htslib' :{'bgzip' : _htslib_config_rule_default.copy()}}

update_config(config_default, config)
config = config_default

rule htslib_bgzip:
    """htslib: run bgzip
    """
    params: cmd = 'bgzip',
            options = config['htslib']['bgzip']['options'],
            runtime = config['htslib']['bgzip']['runtime']
    input: vcf = "{prefix}.vcf"
    output: bgzip = "{prefix}.vcf.gz"
    threads: config['htslib']['bgzip']['threads']
    conda: "env.yaml"
    shell: "{params.cmd} {params.options} -i -@ {threads} {input.vcf}"
