# -*- snakemake -*-
include: 'diamond.settings.smk'

rule diamond:
    """diamond

    Run diamond.
    """
    version: "0.1"
    params: options = config["diamond"]["options"]
    input: 
