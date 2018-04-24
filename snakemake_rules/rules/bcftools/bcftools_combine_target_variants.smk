# -*- snakemake -*-
include: "bcftools.settings.smk"
include: "bcftools_call_targets.smk"
include: "../pybedtools/pybedtools_make_bed_targets.smk"
include: "../gatk/gatk_combine_variants_targets.smk"

rule bcftools_combine_target_variants:
    input: vcf = "{prefix}.bcftools.combined.vcf.gz"
