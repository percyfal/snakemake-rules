# -*- snakemake -*-
include: "utils.settings.smk"
include: "rulegraph.rule"

rule rulegraph_png:
    """Convert rulegraph to png"""
    version: "0.1"
    input: rulegraph = "{prefix}_rulegraph.dot"
    output: png = "{prefix}_rulegraph.png"
    shell: "dot -T png {input.rulegraph} -o {output.png}"

localrules: rulegraph_png
