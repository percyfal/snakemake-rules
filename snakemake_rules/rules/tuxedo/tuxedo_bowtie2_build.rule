# -*- snakemake -*-
include: "tuxedo.settings"

rule tuxedo_bowtie_build2:
    """Bowtie build index"""
    params: ref = config['tuxedo']['ref'],
            cmd = 'bowtie2-build'
    input: "{prefix}" + ".fa"
    output: expand("{{prefix}}{ext}", ext=config['tuxedo']['build']['ext_v2'])
    shell: "{params.cmd} {input} {wildcards.prefix}"
