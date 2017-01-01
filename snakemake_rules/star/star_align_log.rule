# -*- snakemake -*-
include: "star.settings"
include: "star_align_pe.rule"

rule star_align_log:
    input: bam = "{prefix}.Aligned.out.bam"
    output: log = "{prefix}.Log.final.out.bak"
    shell: "cp {wildcards.prefix}.Log.final.out {output.log}"

rule star_align_log_sorted:
    input: bam = "{prefix}.Aligned.sortedByCoord.out.bam"
    output: log = "{prefix}.Log.final.out.bak"
    shell: "cp {wildcards.prefix}.Log.final.out {output.log}"
    

localrules: star_align_log, star_align_log_sorted           
