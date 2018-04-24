# -*- snakemake -*-
include: "utils.settings"
include: "dag.rule"

rule dag_png:
    """Convert dag to png"""
    version: "0.1"
    input: dag = "{prefix}_dag.dot"
    output: png = "{prefix}_dag.png"
    shell: "dot -T png {input.dag} -o {output.png}"
