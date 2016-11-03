# -*- snakemake -*-
include: "tuxedo.settings"

rule tuxedo_bowtie_build:
    """Bowtie build index"""
    params: ref = config['tuxedo']['ref'],
            cmd = 'bowtie-build'
    input: "{prefix}" + ".fa"
    output: expand("{{prefix}}{ext}", ext=config['tuxedo']['build']['ext_v1'])
    shell: "{params.cmd} {input} {wildcards.prefix}"

