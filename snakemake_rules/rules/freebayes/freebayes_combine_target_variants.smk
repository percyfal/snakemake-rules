# -*- snakemake -*-
include: "freebayes.settings.smk"
include: "freebayes_targets.smk"
include: "../htslib/htslib_bgzip.smk"
include: "../pybedtools/pybedtools_make_bed_targets.smk"

if config['freebayes']['freebayes_targets']['make_windows']:
    include: "../gatk/gatk_combine_variants_windows.smk"
else:
    include: "../gatk/gatk_combine_variants_targets.smk"

rule freebayes_combine_target_variants:
    input: vcf = "{prefix}.freebayes.combined.vcf.gz"
