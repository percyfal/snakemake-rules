# -*- snakemake -*-
include: "plink.settings.smk"

rule plink_bed_blocks:
    """calculate linkage disequilibrium blocks(?). 
    """
    params: cmd = config['plink']['cmd'],
            options = config['plink']['options'],
            runtime = config['plink']['runtime']
    input: "{prefix}.bed"
    output: "{prefix}.blocks"
    threads: config['plink']['threads']
    shell: "{params.cmd} {params.options}  --bfile {wildcards.prefix} --blocks --out {wildcards.prefix}"

