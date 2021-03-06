# -*- snakemake -*-
include: "rseqc.settings.smk"

RULES = ['rseqc_clipping_profile',
         'rseqc_geneBody_coverage',
         'rseqc_junction_annotation',
         'rseqc_read_GC',
         'rseqc_read_NVC',
         'rseqc_read_distribution',
         'rseqc_read_duplication',
         'rseqc_read_quality']

for rule in RULES:
    if not rule in workflow._rules:
        include: rule + ".smk"

rule rseqc_qc_8:
    """Run 8 RSeQC commands on a bam file"""
    input: read_dup = "{prefix}_rseqc/read_dup.pos.DupRate.xls", 
           clipping_profile = "{prefix}_rseqc/clippingprofile.clipping_profile.xls",
           geneBody_coverage = "{prefix}_rseqc/geneBody_coverage.geneBodyCoverage.txt",
           junction_annotation = "{prefix}_rseqc/junction_annotation_refseq.txt",
           read_GC = "{prefix}_rseqc/read_GC.GC.xls",
           read_NVC = "{prefix}_rseqc/read_NVC.NVC.xls",
           read_quality = "{prefix}_rseqc/read_quality.qual.r",
           read_distribution = "{prefix}_rseqc/read_distribution.txt"
    output: "{prefix}_rseqc/rseqc_qc_8.txt"
    shell: "echo `date` > {output}"


localrules: rseqc_qc_8           
