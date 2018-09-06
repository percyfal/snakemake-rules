# -*- snakemake -*-
include: "freebayes.settings"
include: "../htslib/htslib_bgzip.rule"

config_default = {
    'freebayes' : {
        'options' : '',
    },
}

update_config(config_default, config)
config = config_default

rule freebayes:
    """Run freebayes on a set of bam files"""
    wildcard_constraints:
        suffix = "(.vcf.gz|.vcf)"
    params: cmd = config['freebayes']['cmd'],
            options = config['freebayes']['options'],
            runtime = config['freebayes']['runtime'],
            bgzip = 'bgzip'
    input: fofn = "{prefix}.bam.fofn", ref = config['freebayes']['ref']
    output: vcf = "{prefix}.freebayes{suffix}"
    threads: config['freebayes']['threads']
    conda: "env.yaml"
    shell:
        "if [[ '{wildcards.suffix}' == 'vcf.gz' ]]; then \n"
        "  {params.cmd} {params.options} -f {input.ref} -L {input.fofn} | {params.bgzip} -c > {output.vcf}\n"
        "  else\n"
        "  {params.cmd} {params.options} -f {input.ref} -L {input.fofn} > {output.vcf}\n"
        "fi"
