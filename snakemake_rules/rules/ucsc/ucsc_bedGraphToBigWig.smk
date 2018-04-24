# -*- snakemake -*-
include: "ucsc.settings.smk"
include: "ucsc_fetchChromSizes.smk"

config_default = {'ucsc': {'bedGraphToBigWig' : _ucsc_config_rule_default.copy()}}
config_default['ucsc']['bedGraphToBigWig'].update({'cmd' : 'bedGraphToBigWig'})

update_config(config_default, config)
config = config_default


rule ucsc_bedGraphToBigWig:
    """Convert bedGraph file to bigWig.

    Run bedGraphToBigWig to convert bedGraph file to bigWig.
    """
    params: cmd = config['ucsc']['bedGraphToBigWig']['cmd'],
            options = config['ucsc']['bedGraphToBigWig']['options'],
            runtime = config['ucsc']['bedGraphToBigWig']['runtime']
    input: wig = "{prefix}.bedGraph",
           sizes = "chrom.sizes"
    output: bigwig = "{prefix}.bigWig"
    threads: config['ucsc']['bedGraphToBigWig']['threads']
    conda: "envs/ucsc_bedgraphtobigwig.yaml"
    shell: "{params.cmd} {params.options} {input.wig} {input.sizes} {output.bigwig}"
