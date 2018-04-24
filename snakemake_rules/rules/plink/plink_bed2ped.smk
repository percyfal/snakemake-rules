# -*- snakemake -*-
include: "plink.settings.smk"

rule plink_bed2ped:
    """convert ped file to binary ped. 
    """
    params: cmd = config['plink']['cmd'],
            options = config['plink']['options'],
            runtime = config['plink']['runtime']
    input: "{prefix}.ped"
    output: "{prefix}.bed"
    threads: config['plink']['threads']
    shell: "{params.cmd} {params.options}  --file {wildcards.prefix} --make-bed --out {wildcards.prefix}"

