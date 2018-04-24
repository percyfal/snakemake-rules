# -*- snakemake -*-
include: "tuxedo.settings.smk"

rule tuxedo_cufflinks_quant:
    """Run cufflinks quantification"""
    params: cmd = config['tuxedo']['cufflinks']['cmd'],
            options = config['tuxedo']['cufflinks']['options'],
            hits = os.path.join("{prefix}.tophat2" if config['tuxedo']['version2'] else "{prefix}.tophat", "accepted_hits.bam"),
    input: tophat = "{prefix}.tophat2" if config['tuxedo']['version2'] else "{prefix}.tophat",
           annot_gtf = config['tuxedo']['cufflinks']['transcript_annot_gtf']
    output: quant = "{prefix}.cufflinks_quant"
    shell: "{params.cmd} {params.options} --GTF {input.annot_gtf} {params.hits} -o {output.quant} &> {output.quant}.log"

