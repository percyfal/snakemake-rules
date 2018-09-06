# -*- snakemake -*-
include: "picard.settings.smk"

rule picard_create_sequence_dictionary_awk:
    """Picard: create interval list dict using awk"""
    params: cmd = config['picard']['cmd'] + CREATE_SEQUENCE_DICTIONARY,
            options = config['picard']['options']
    input: bed="{prefix}.bed", dict=config['picard']['ref'].replace(".fa", ".dict")
    output: "{prefix}.dict"
    conda: "env.yaml"
    shell: config['comp']['cat'] + " {input.dict} > {output}; " + config['comp']['awk'] + " '{{printf(\"%s\\t%s\\t%s\\t%s\\t%s\\n\", $1,$2,$3,\"+\",$4); FS=\"\t\"}}' {input.bed} >> {output}"

