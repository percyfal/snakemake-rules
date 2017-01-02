# -*- snakemake -*-
include: "emacs.settings"

config_default = {
    'emacs' : {
        'org_to_reveal' : {
            'options' : '',
            'config' : 'config.el',
        },
    },
}

update_config(config_default, config)
config = config_default


rule emacs_org_to_reveal:
    """Convert org file to reveal presentation"""
    version: "0.1"
    params: options = config['emacs']['org_to_reveal']['options'],
            cmd = config['emacs']['cmd']
    input: org = "{prefix}.org", config = config['emacs']['org_to_reveal']['config']
    output: html = "{prefix}.html"
    shell: "{params.cmd} --batch --no-init-file --load $(basename {input.config}) --find-file {input.org} --funcall org-reveal-export-to-html"

localrules: org_to_reveal
