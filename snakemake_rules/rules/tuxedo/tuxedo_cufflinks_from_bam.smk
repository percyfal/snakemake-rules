# -*- snakemake -*-
include: "tuxedo.settings.smk"

rule tuxedo_cufflinks_from_bam:
    """Run cufflinks on bam file"""
    params: cmd = config['tuxedo']['cufflinks']['cmd'],
            options = config['tuxedo']['cufflinks']['options']
    input: bam="{prefix}.bam"
    output: "{prefix}.cufflinks"
    shell: "{params.cmd} {params.options} {input.bam} -o {output}.tmp &> {output}.log && mv {output}.tmp {output} && mv {output}.log {output}"
