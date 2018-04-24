# -*- snakemake -*-
include: "plink.settings.smk"

rule plink_chr_bed:
    """generate chromosome-wise bed files. 
    """
    params: cmd = config['plink']['cmd'],
            options = config['plink']['options'],
            chr = config['plink']['chr'],
            runtime = config['plink']['runtime'],
    input: "{prefix}.bed"
    output: "{prefix}.chr{params.chr}.bed"
    threads: config['plink']['threads']
    shell: "{params.cmd} {params.options} --bfile {wildcards.prefix} --make-bed --out {wildcards.prefix}.chr{params.chr}"
