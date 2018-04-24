include: 'r.settings'

config_default = {'r' :{'rmarkdown_bookdown' : _r_config_rule_default.copy()}}
config_default['r']['rmarkdown_bookdown']['gitbook'] = {'options': "split_by=\"chapter\""}
config_default['r']['rmarkdown_bookdown']['tufte_html_book'] = {'options': "split_by=\"chapter\""}
config_default['r']['rmarkdown_bookdown']['html_book'] = {'options': "split_by=\"chapter\""}
config_default['r']['rmarkdown_bookdown']['tufte_html2'] = {'options': ""}
config_default['r']['rmarkdown_bookdown']['pdf_book'] = {'options': ""}

update_config(config_default, config)
config = config_default

rule rmarkdown_bookdown:
    """Convert Rmarkdown to bookdown output based on file suffix"""
    wildcard_constraints:
        suffix = "(pdf|html)",
        booktype = "(gitbook|tufte_html_book|html_book|tufte_html2|pdf_book|pdf_document2)"
    params: runtime = config['r']['rmarkdown_bookdown']['runtime'],
            xvfb = "" if config['r']['X11'] else "xvfb-run -a ",
            options = ",{}".format(config['r']['rmarkdown_bookdown']['options']) if config['r']['rmarkdown_bookdown']['options'] else "",
            cmd = config['r']['script']
    input: rmd = "{prefix}.Rmd",
           result_files = config['__RESULT_FILES__']
    output: out = "{prefix}.bookdown.{booktype}.{suffix}"
    threads: config['r']['rmarkdown_bookdown']['threads']
    run:
        wd = os.path.dirname(os.path.abspath(input.rmd))
        infile = os.path.basename(input.rmd)
        outfile = os.path.basename(output.out)
        bookdown_format = wildcards.booktype
        bookdown_options = config['r']['rmarkdown_bookdown'][bookdown_format]['options']

        command = "{xvfb} {cmd} --verbose -e 'bookdown::render_book(\"{infile}\", output_format=bookdown::{format}({bookdown_options}), output_file=\"{out}\"{options})'".format(
                xvfb=params.xvfb, cmd=params.cmd, infile=infile, out=outfile, options=params.options, format=bookdown_format, bookdown_options=bookdown_options)
        with cd(wd, logger):
            if os.path.exists("_main.Rmd"):
                os.unlink("_main.Rmd")
            shell(command)
            ## Assume index or abstract are main pages for now
            indices = ["_book/index.html", "_book/abstract.html"]
            for fn in indices:
                if os.path.exists(fn) and not os.path.exists(os.path.basename(output.out)):
                    print("linking output {} to index file {}".format(output.out, fn))
                    shell("ln -s {} {}".format(fn, os.path.basename(output.out)))
                    break
