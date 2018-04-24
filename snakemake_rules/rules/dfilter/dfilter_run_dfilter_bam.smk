# -*- snakemake -*-
include: "dfilter.settings.smk"

rule dfilter_run_dfilter_bam:
    """Run run_dfilter command. Currently only works on one file."""
    params: options = config['dfilter']['options'],
            cmd = config['dfilter']['cmd']
    input: chipfile = "{prefix}.bam"
    output: bed = "{prefix}.dfilt.bed",
            wig = "{prefix}.dfilt.bed.wig" if "-wig" in config['dfilter']['options'] else []
    shell: "{params.cmd} {params.options} -f=bam -d={input.chipfile} -o={output.bed}"

