# -*- snakemake -*-
include: "utils.settings"

rule gzip:
    """gzip"""
    version: "0.1"
    input: infile = "{prefix}"
    output: outfile = "{prefix}.gz"
    shell: "gzip -v {input.infile}"
