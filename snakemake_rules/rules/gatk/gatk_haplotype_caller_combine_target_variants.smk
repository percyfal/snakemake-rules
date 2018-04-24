# -*- snakemake -*-
include: "gatk.settings.smk"
include: "gatk_haplotype_caller_targets.rule"
include: "../pybedtools/pybedtools_make_bed_targets.rule"
include: "../gatk/gatk_combine_gvcfs_targets.rule"
include: "../gatk/gatk_combine_variants_targets.rule"

rule gatk_haplotype_caller_combine_target_variants:
    wildcard_constraints:
        mode = "(.g|)",
        suffix = "(.vcf.gz|.vcf)"
    input: vcf = "{prefix}.haplotype_caller.combined{mode}{suffix}",
           tbi = "{prefix}.haplotype_caller.combined{mode}{suffix}.tbi"
