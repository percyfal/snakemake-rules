# -*- snakemake -*-
include: "jellyfish.settings.smk"
include: "jellyfish_count.smk"

config_default = {
    'jellyfish': {
        'histo': {
            'options': "",
            'runtime' : config['jellyfish']['runtime'],
        },
    },
}

update_config(config_default, config)
config = config_default

rule jellyfish_histo:
    params: cmd = config['jellyfish']['cmd'],
            options = config['jellyfish']['histo']['options'],
            runtime = config['jellyfish']['histo']['runtime'],
    input: counts = "{prefix}.mer_counts.jf"
    output: histo = "{prefix}.histo"
    threads: config['jellyfish']['threads']
    conda: "env.yaml"
    shell:
        "{params.cmd} histo {params.options} {input.counts} -o {output.histo}"
