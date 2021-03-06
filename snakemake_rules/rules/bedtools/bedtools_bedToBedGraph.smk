# -*- snakemake -*-
include: "bedtools.settings.smk"
include: join(os.pardir, "ucsc", "ucsc_fetchChromSizes.smk")

config_default = {'bedtools' :{'genomecov' : _bedtools_config_rule_default.copy()}}

update_config(config_default, config)
config = config_default


rule bedtools_bedToBedGraph:
    """Convert bed to bedGraph
    """
    params: cmd = " ".join([BEDTOOLS, BEDTOOLS_GENOMECOV]),
            options=config['bedtools']['genomecov']['options'],
            runtime=config['bedtools']['genomecov']['runtime'],
    input: infile = "{prefix}.bed",
           chrom = config['bedtools']['chromsizes']
    output: bedgraph = "{prefix}.bedGraph"
    threads: config['bedtools']['genomecov']['threads']
    conda: "env.yaml"
    shell: "{params.cmd} {params.options} -bg -trackline -trackopts name={wildcards.prefix} -i {input.infile} -g {input.chrom} > {output.bedgraph}"
