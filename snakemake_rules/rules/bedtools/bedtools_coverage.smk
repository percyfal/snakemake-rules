# -*- snakemake -*-
include: "bedtools.settings.smk"

config_default = {'bedtools' :{'coverage' : _bedtools_config_rule_default.copy()}}

update_config(config_default, config)
config = config_default

rule bedtools_coverage:
    """Returns the depth and breadth of coverage of features from B
    on the intervals in A.
    """
    wildcard_constraints:
        suffix = "(bam|fofn)",
        window = "\d*"
    params: cmd = " ".join([BEDTOOLS, BEDTOOLS_COVERAGE]),
            windowcmd = " ".join([BEDTOOLS, BEDTOOLS_MAKEWINDOWS]),
            options=config['bedtools']['coverage']['options'],
            runtime = config['bedtools']['coverage']['runtime'],
    input: a = config['bedtools']['coverage']['afile'], b = "{prefix}.{suffix}"
    output: bed = "{prefix}.{suffix}.coverage{window}.bed"
    threads: config['bedtools']['coverage']['threads']
    conda: "env.yaml"
    shell:
        "if [ \"{wildcards.window}\" = \"\" ]; then\n"
        "    windowcmd=\"cat {input.a}\"\n"
        "else\n"
        "    windowcmd=\"{params.windowcmd} -b {input.a} -w {wildcards.window}\"\n"
        "fi\n"
        "bedfiles=\"\"; if [ \"{wildcards.suffix}\" = \"bam\" ]; then\n"
        "   bedfiles=\"-b {input.b}\"\n"
        "else\n"
        "   bedfiles=\"$(cat {input.b} | sed -e 's/[^ ][^ ]*/-b &/g' | tr '\\n' ' ')\";\n"
        "fi\n"
        "${{windowcmd}} |  {params.cmd} {params.options} -a stdin ${{bedfiles}} > {output}"

