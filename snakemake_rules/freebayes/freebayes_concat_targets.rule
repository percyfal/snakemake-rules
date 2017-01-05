# -*- snakemake -*-
include: "freebayes.settings"
include: "freebayes_targets.rule"
include: "../bcftools/bcftools.settings"
include: "../pybedtools/pybedtools_make_bed_targets.rule"
include: "../htslib/htslib_bgzip.rule"
include: "../samtools/samtools_tabix_vcf.rule"

config_default = {
    'freebayes' : {
        'concat_targets': {
            'options' : "",
            'runtime' : "01:00:00",
        },
    },
}

update_config(config_default, config)
config = config_default

rule freebayes_concat_targets:
    """Concatenate freebayes vcfs"""
    params: cmd = config['bcftools']['cmd'],
            options = config['freebayes']['concat_targets']['options'],
            runtime = config['freebayes']['concat_targets']['runtime']
    input: vcf = ["{{prefix}}.freebayes.{partition}.vcf.gz".format(partition=p+1) for p in range(config['pybedtools']['make_bed_targets']['partitions'])]
    output: vcf = "{prefix}.freebayes.concat.vcf.gz"
    threads: config['freebayes']['threads']
    shell: "{params.cmd}  concat {params.options} -O z {input.vcf} -o {output.vcf}"


