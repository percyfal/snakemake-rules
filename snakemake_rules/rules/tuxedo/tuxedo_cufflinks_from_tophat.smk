# -*- snakemake -*-
include: "tuxedo.settings"

rule tuxedo_cufflinks_from_tophat:
    """Run cufflinks based on tophat output directory"""
    params: cmd = config['tuxedo']['cufflinks']['cmd'],
            options = config['tuxedo']['cufflinks']['options']
    input: tophat = os.path.join("{prefix}.tophat2", "accepted_hits.bam") if config['tuxedo']['version2'] else os.path.join("{prefix}.tophat", "accepted_hits.bam")
    output: "{prefix}.cufflinks"
    shell: "{params.cmd} {params.options} {input.tophat} -o {output}.tmp &> {output}.log && mv {output}.tmp {output} && mv {output}.log {output}"

