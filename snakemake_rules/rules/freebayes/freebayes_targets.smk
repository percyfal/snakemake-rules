# -*- snakemake -*-
include: "freebayes.settings.smk"
include: "../samtools/samtools_faidx.rule"
include: "../htslib/htslib_bgzip.rule"

config_default = {'freebayes' : {'freebayes_targets' : _freebayes_config_rule_default.copy()}}
config_default['freebayes']['freebayes_targets']['options'] = config['freebayes']['options']
config_default['freebayes']['freebayes_targets']['make_windows'] = False

update_config(config_default, config)
config = config_default

if config['freebayes']['freebayes_targets']['make_windows']:
    include: "../pybedtools/pybedtools_make_bed_windows.rule"
    npartitions = config['pybedtools']['make_bed_windows']['partitions']
else:
    include: "../pybedtools/pybedtools_make_bed_targets.rule"
    npartitions = config['pybedtools']['make_bed_targets']['partitions']

rule freebayes_targets:
    """Run freebayes on a set of bam files on targets given in bed file.

    Takes as input a fofn containing a list of bam file names, a
    reference, the reference index, and a bed file containing target
    regions in which variant calling is performed. The targets
    wildcard is often a number, resulting from the partitioning of a
    larger bed files into several bed files of subregions.

    """
    wildcard_constraints:
        window = "(|window-)",
        targets = "\d+",
        suffix = "(.vcf.gz|.vcf)"
    params: cmd = config['freebayes']['cmd'],
            options = config['freebayes']['freebayes_targets']['options'],
            runtime = config['freebayes']['freebayes_targets']['runtime'],
            bgzip = 'bgzip'
    input: fofn = "{prefix}.fofn",
           ref = config['freebayes']['freebayes_targets']['ref'],
           fai = config['freebayes']['freebayes_targets']['ref'] + ".fai",
           targets = "{prefix}.{window}{targets}.bed"
    output: vcf = "{prefix}.freebayes.{window}{targets}{suffix}"
    threads: config['freebayes']['freebayes_targets']['threads']
    shell:
        "if [[ '{wildcards.suffix}' == '.vcf.gz' ]]; then \n"
        "  {params.cmd} {params.options} -f {input.ref} -L {input.fofn} -t {input.targets} | {params.bgzip} > {output.vcf}\n"
        "else\n"
        "  {params.cmd} {params.options} -f {input.ref} -L {input.fofn} -t {input.targets} > {output.vcf}\n"
        "fi"

ruleorder: freebayes_targets > htslib_bgzip
