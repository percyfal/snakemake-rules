# -*- snakemake -*-
include: "plink.settings.smk"

rule plink_ped2ld_r2:
    """calculate r-squared linked disequilibrium from ped file. 
    """
    params: cmd = config['plink']['cmd'],
            options = config['plink']['options'],
            runtime = config['plink']['runtime']
    input: "{prefix}.ped"
    output: "{prefix}.r2.ld"
    threads: config['plink']['threads']
    shell: "{params.cmd} {params.options}  --file {wildcards.prefix} --r2 --out {wildcards.prefix}.r2"

