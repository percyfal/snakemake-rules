# -*- snakemake -*-
include: "freebayes.settings.smk"
include: "freebayes_targets.rule"
include: "../htslib/htslib_bgzip.rule"
include: "../pybedtools/pybedtools_make_bed_targets.rule"

if config['freebayes']['freebayes_targets']['make_windows']:
    include: "../gatk/gatk_combine_variants_windows.rule"
else:
    include: "../gatk/gatk_combine_variants_targets.rule"

rule freebayes_combine_target_variants:
    input: vcf = "{prefix}.freebayes.combined.vcf.gz"
