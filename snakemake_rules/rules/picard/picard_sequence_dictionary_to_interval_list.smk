# -*- snakemake -*-
include: "picard.settings.smk"

rule picard_sequence_dictionary_to_interval_list:
    """Picard: create interval list dict from sequence dictionary"""
    params: runtime = "01:00:00"
    input: dict = "{prefix}.dict"
    output: interval_list = "{prefix}.interval_list"
    conda: "env.yaml"
    threads: 1
    shell:
        "cat {input.dict} > {output.interval_list};"
        "cat {input.dict} | sed -e \"s/\(LN:\|SN:\)//g\" | grep \"@SQ\" | awk '{{printf(\"%s\\t%s\\t%s\\t%s\\t%s\\n\", $2, 1,$3,\"+\",\".\"); FS=\"\t\"}}' >> {output.interval_list}"

localrules: picard_sequence_dictionary_to_interval_list
