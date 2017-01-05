# -*- snakemake -*-
# See https://groups.google.com/forum/#!topic/rna-star/tvajn49WTYk for
# setting up RSEM with STAR alignments
include: "rsem.settings"

rule rsem_calculate_expression:
    """Calculate RSEM expression from bam"""
    params: cmd = config['rsem']['calculate-expression']['cmd'],
            options = " ".join(["--bam",
                                config['rsem']['calculate-expression']['options']]),
            index = str(config['rsem']['index']),
            runtime = config['rsem']['runtime']
    input: index = str(config['rsem']['index']) + config['rsem']['ref_sfx'],
           bam = "{prefix}.bam"
    output: isoforms = "{prefix}.isoforms.results", genes = "{prefix}.genes.results"
    threads: config['rsem']["threads"]
    conda: "env.yaml"
    shell: "{params.cmd} {params.options} -p {threads} {input.bam} {params.index} {wildcards.prefix}"
