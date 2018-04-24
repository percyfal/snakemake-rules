# -*- snakemake -*-
include: 'r.settings.smk'

config_default = {'r' :{'rmarkdown_revealjs' : _r_config_rule_default.copy()}}
config_default['r']['rmarkdown_revealjs']['css'] = False

update_config(config_default, config)
config = config_default

rule rmarkdown_revealjs:
    """Create revealjs presentation from Rmarkdown"""
    params: runtime = config['r']['rmarkdown_revealjs']['runtime'],
            xvfb = "" if config['r']['X11'] else "xvfb-run -a ",
            cmd = config['r']['script']
    input: rmd = "{prefix}.Rmd",
           css = "{prefix}.css" if config_default['r']['rmarkdown_revealjs']['css'] else []
    output: html = "{prefix}.revealjs{standalone, (|.standalone)}.html"
    threads: config['r']['rmarkdown_revealjs']['threads']
    run:
        wd = os.path.dirname(os.path.abspath(input.rmd))
        infile = os.path.basename(input.rmd)
        outfile = os.path.basename(output.html)
        standalone = "true"  if wildcards.standalone == ".standalone" else "false"
        command = "{xvfb} {cmd} -e 'rmarkdown::render(\"{infile}\", output_file=\"{out}\", output_options=(standalone=\"{standalone}\"))'".format(
            xvfb=params.xvfb, cmd=params.cmd, infile=infile, out=outfile, standalone=standalone)
        with cd(wd, logger):
            if os.path.exists("_main.Rmd"):
                os.unlink("_main.Rmd")
            shell(command)

    # shell:
    #     "standalone=\"false\"; if [  \"{wildcards.standalone}\" = \".standalone\" ]; then standalone=\"true\"; fi\n"
    #     "{params.xvfb} {params.cmd} -e 'rmarkdown::render(\"{wildcards.prefix}{wildcards.standalone}.Rmd\", output_file=\"{wildcards.prefix}.revealjs{wildcards.standalone}.html\", output_options=(standalone=\"${{standalone}}\"))'"
