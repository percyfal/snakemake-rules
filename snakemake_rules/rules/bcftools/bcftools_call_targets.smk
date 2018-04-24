# -*- snakemake -*-
include: "bcftools.settings.smk"
include: "../pybedtools/pybedtools_make_bed_targets.smk"

config_default = {'bcftools' :{'call_targets' : _bcftools_config_rule_default.copy()}}
config_default['bcftools']['call_targets'].update(
    {
        'options' : '-vm',
        'runtime' : '100:00:00',
        'mpileup_options' : "-ug",
    })


update_config(config_default, config)
config = config_default

rule bcftools_call_targets:
    """bcftools call targets: call variants from samtools mpileup on target region

    Takes as input a fofn containing a list of bam file names, a
    reference, the reference index, and a bed file containing target
    regions in which variant calling is performed. The targets
    wildcard is often a number, resulting from the partitioning of a
    larger bed files into several bed files of subregions.

    """
    params: cmd = config['bcftools']['cmd'],
            options = config['bcftools']['call_targets']['options'],
            mpileup_options = config['bcftools']['call_targets']['mpileup_options'],
            runtime = config['bcftools']['call_targets']['runtime']
    input: fofn = "{prefix}.fofn",
           ref = config['bcftools']['ref'],
           fai = config['bcftools']['ref'] + ".fai",
           targets = "{prefix}.{targets}.bed"
    output: vcf = "{prefix}.bcftools.{targets}.vcf.gz"
    threads: config['bcftools']['call_targets']['threads']
    conda: "env.yaml"
    shell: "samtools mpileup {params.mpileup_options} -f {input.ref} -b {input.fofn} -l {input.targets} | {params.cmd} call {params.options} -O z -o {output.vcf}"

