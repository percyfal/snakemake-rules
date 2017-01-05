# -*- snakemake -*-
include: "homer.settings"

config_default = {'homer' :{'annotatePeaks' : _homer_config_rule_default.copy()}}
config_default['homer']['annotatePeaks'].update({'gtf': None})

update_config(config_default, config)
config = config_default

options =  " ".join([
    config['homer']['annotatePeaks']['options'],
    "-gtf {}".format(config['homer']['annotatePeaks']['gtf']) if config['homer']['annotatePeaks']['gtf'] else ""])

gtf = [config['homer']['annotatePeaks']['gtf']] if config['homer']['annotatePeaks']['gtf'] else []

rule annotatePeaks:
    """homer: run annotatePeaks.pl for annotating peaks."""
    params: cmd = "annotatePeaks.pl",
            build = config['homer']['genome_build'],
            options = options,
            runtime = config['homer']['annotatePeaks']['runtime'],
    input: bed = "{prefix}.bed",
           gtf = gtf
    log: "{prefix}.annotatePeaks.log"
    output: tsv = "{prefix}.annotatePeaks.tsv",
            stats = "{prefix}.annotatePeaks.stats.tsv"
    threads: config['homer']['annotatePeaks']['threads']
    shell: "{params.cmd} {input.bed} {params.build} -annStats {output.stats} {params.options} > {output.tsv} 2> {log}"
