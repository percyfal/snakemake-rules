# -*- snakemake -*-
include: "bwa.settings.smk"

rule bwa_link_ref:
    """bwa link reference file to bwa index directory"""
    input: ref = config['bwa']['ref']
    output: indexref = config['bwa']['index']
    shell: "if [ ! -e {output.indexref} ]; then ln -s {input.ref} {output.indexref}; fi"


localrules: bwa_link_ref           
