# -*- snakemake -*-
include: 'diamond.settings.smk'

# See http://megan.informatik.uni-tuebingen.de/t/generic-pipeline-using-diamond-and-megan6/50
# Adding output paths here as that seems to be common practice
config_default = {
    'diamond' : {
        'blastx' : {
            'options' : "",
            'fastq' : "",
            'daa' : "",
        }
    }
}

update_config(config_default, config)
config = config_default


rule diamond_blastx:
    """diamond_blastx

    Run diamond blastx.
    """
    version: "0.1"
    params: options = config["diamond"]["options"]
    input: join(config["diamond"]["blastx"]["fastq"], "{prefix}")
