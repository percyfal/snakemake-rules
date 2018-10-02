# -*- snakemake -*-
include: "freebayes.settings.smk"
include: "freebayes_targets.smk"
include: "../bcftools/bcftools_index.smk"
include: "../pybedtools/pybedtools_make_bed_targets.smk"
include: "../htslib/htslib_bgzip.smk"

config_default = {
    'freebayes' : {
        'merge_targets': {
            'options' : "",
            'runtime' : "01:00:00",
        },
    },
}

update_config(config_default, config)
config = config_default

rule freebayes_merge_targets:
    """Merge freebayes vcfs"""
    params: cmd = config['bcftools']['cmd'],
            options = config['freebayes']['merge_targets']['options'],
            runtime = config['freebayes']['merge_targets']['runtime']
    input: vcf = ["{{prefix}}.freebayes.{partition}.vcf.gz".format(partition=p+1) for p in range(config['pybedtools']['make_bed_targets']['partitions'])],
           csi = ["{{prefix}}.freebayes.{partition}.vcf.gz.csi".format(partition=p+1) for p in range(config['pybedtools']['make_bed_targets']['partitions'])]
    output: vcf = "{prefix}.freebayes.merge.vcf.gz"
    threads: config['freebayes']['threads']
    shell: "{params.cmd}  merge {params.options} -O z {input.vcf} -o {output.vcf}"


