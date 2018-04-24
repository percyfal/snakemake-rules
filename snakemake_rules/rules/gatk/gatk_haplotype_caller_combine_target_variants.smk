# -*- snakemake -*-
include: "gatk.settings.smk"
include: "gatk_haplotype_caller_targets.smk"
include: "../pybedtools/pybedtools_make_bed_targets.smk"
include: "../gatk/gatk_combine_gvcfs_targets.smk"
include: "../gatk/gatk_combine_variants_targets.smk"

rule gatk_haplotype_caller_combine_target_variants:
    wildcard_constraints:
        mode = "(.g|)",
        suffix = "(.vcf.gz|.vcf)"
    input: vcf = "{prefix}.haplotype_caller.combined{mode}{suffix}",
           tbi = "{prefix}.haplotype_caller.combined{mode}{suffix}.tbi"
