# -*- snakemake -*-
include: "bcftools.settings.smk"

config_default = {'bcftools' :{'call' : _bcftools_config_rule_default.copy()}}
config_default['bcftools']['call'].update(
    {
        'options' : '-vmO z',
        'runtime' : '100:00:00',
        'mpileup_options' : "-ug",
    })

update_config(config_default, config)
config = config_default

rule bcftools_call:
    """bcftools call: call variants from samtools mpileup"""
    params: cmd = config['bcftools']['cmd'],
            options = config['bcftools']['call']['options'],
            mpileup_options = config['bcftools']['call']['mpileup_options'],
            runtime = config['bcftools']['call']['runtime']
    input: fofn = "{prefix}.bam.fofn", ref = config['bcftools']['ref']
    output: vcf = "{prefix}.bcftools.vcf.gz"
    threads: config['bcftools']['call']['threads']
    conda: "env.yaml"
    shell: "samtools mpileup {params.mpileup_options} -f {input.ref} -b {input.fofn} | {params.cmd} call {params.options} -o {output.vcf}"

