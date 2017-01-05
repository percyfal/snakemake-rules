# -*- snakemake -*-
include: "htslib.settings"

config_default = {'htslib' :{'plot-vcfcheck' : _htslib_config_rule_default.copy()}}

update_config(config_default, config)
config = config_default

rule htslib_plot_vcfcheck:
    """htslib: plot vcf check outut file
    
    Output from vcftools stats.
    """
    params: cmd='plot-vcfcheck',
            options = config['htslib']['plot-vcfcheck']['options'],
            runtime = config['htslib']['plot-vcfcheck']['runtime']
    input: "{prefix}.chk"
    output: "{prefix}-summary.pdf"
    threads: config['htslib']['plot-vcfcheck']['threads']
    conda: "env.yaml"
    shell: "{params.cmd} {input} -p {wildcards.prefix}"
