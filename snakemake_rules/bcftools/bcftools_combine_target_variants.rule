# -*- snakemake -*-
include: "bcftools.settings"
include: "bcftools_call_targets.rule"
include: "../pybedtools/pybedtools_make_bed_targets.rule"
include: "../gatk/gatk_combine_variants_targets.rule"

rule bcftools_combine_target_variants:
    input: vcf = "{prefix}.bcftools.combined.vcf.gz"
