# -*- snakemake -*-
include: "ucsc.settings"

rule ucsc_genepred_to_refFlat:
    """Convert genepred to refFlat.

    See
    https://github.com/chapmanb/cloudbiolinux/blob/master/utils/prepare_tx_gff.py#L355;
    we need to add extra column at start.

    """
    params: runtime = config['ucsc']['runtime']
    input: genepred = "{prefix}.genePred"
    output: refflat = "{prefix}.refFlat"
    threads: 1
    shell: "awk '{{printf(\"%s\\t%s\\n\", $1, $0)}}'  {input.genepred} > {output.refflat}"

localrules: ucsc_genepred_to_refFlat           
