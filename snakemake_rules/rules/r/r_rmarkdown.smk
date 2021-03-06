# -*- snakemake -*-
include: 'r.settings.smk'

config_default = {'r' :{'rmarkdown' : _r_config_rule_default.copy()}}

update_config(config_default, config)
config = config_default

rule rmarkdown:
    """Convert Rmarkdown to output based on file suffix"""
    wildcard_constraints:
        suffix = "(pdf|html)"
    params: runtime = config['r']['rmarkdown']['runtime'],
            xvfb = "" if config['r']['X11'] else "xvfb-run -a ",
            options = ",{}".format(config['r']['rmarkdown']['options']) if config['r']['rmarkdown']['options'] else "",
            cmd = config['r']['script']
    input: rmd = "{prefix}.Rmd"
    output: out = "{prefix}.{suffix}"
    threads: config['r']['rmarkdown']['threads']
    shell:
        "{params.xvfb} {params.cmd} -e 'rmarkdown::render(\"{wildcards.prefix}.Rmd\", output_format=\"{wildcards.suffix}_document\", output_file=\"{output.out}\"{params.options})'"
