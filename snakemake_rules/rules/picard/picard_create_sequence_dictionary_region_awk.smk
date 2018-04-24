# -*- snakemake -*-
include: "picard.settings.smk"

rule picard_create_region_dict_awk:
    """Picard: create interval list dict for region using awk"""
    params: cmd = config['picard']['cmd'] + CREATE_SEQUENCE_DICTIONARY,
            options = config['picard']['options']
    input: bed="{prefix}.region_{gene}.bed", dict=config['picard']['ref'].replace(".fa", ".dict")
    output: "{prefix}.region_{gene}.dict"
    conda: "env.yaml"
    shell: config['comp']['cat'] + " {input.dict} > {output}; " + config['comp']['awk'] + " '{{printf(\"%s\\t%s\\t%s\\t%s\\t%s\\n\", $1,$2,$3,\"+\",$4); FS=\"\t\"}}' {input.bed} >> {output}"
