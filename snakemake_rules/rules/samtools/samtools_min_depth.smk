# -*- snakemake -*-
include: "samtools.settings.smk"
include: "../bedtools/bedtools.settings.smk"

config_default = {'samtools' :{'min_depth' : _samtools_config_rule_default.copy()}}

update_config(config_default, config)
config = config_default


rule samtools_min_depth:
    """Run samtools depth, exporting only regions with a given depth

    Outputs a region-merged bed file.
    """
    params: cmd = config['samtools']['cmd'],
            bedtoolscmd = " ".join([BEDTOOLS, BEDTOOLS_MERGE]),
            options = config['samtools']['min_depth']['options'],
            runtime = config['samtools']['min_depth']['runtime']
    input: fofn = "{prefix}.fofn"
    output: bed = "{prefix}.above{mincov}.bed"
    threads: config['samtools']['min_depth']['threads']
    conda: "env.yaml"
    shell: "{params.cmd} {params.options} depth -f {input.fofn} "
           " | awk 'BEGIN {{OFS=\"\\t\"}} {{ {{a=0; for (i=3; i<=NF; i++) a+=$i}}; if (a>={wildcards.mincov}) {{print $1, $2 - 1, $2, a;}} }}' "
           " | {params.bedtoolscmd} -i - > {output.bed}"
