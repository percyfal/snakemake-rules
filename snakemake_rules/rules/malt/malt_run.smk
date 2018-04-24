# -*- snakemake -*-
include: 'malt.settings.smk'

config_default = {'malt' :{'run' : _malt_config_rule_default.copy()}}
config_default['malt']['run']['mode'] = "BlastX"

update_config(config_default, config)
config = config_default

rule malt_run:
    """Align sequences using MALT (MEGAN alignment tool)

    The input is a fasta or fastq file. Note that MALT does not have a
    paired-end mode. See
    http://megan.informatik.uni-tuebingen.de/t/paired-end-reads-in-malt-and-megan/391.

    """
    wildcard_constraints:
        fa = "(.fasta.gz|.fastq.gz)"
    params: runtime = config['malt']['run']['runtime'],
            xvfb = "" if config['malt']['X11'] else "xvfb-run -a ",
            options = config['malt']['run']['options'],
            mode = config['malt']['run']['mode']
    input: fa = "{prefix}{ext}", index = config['malt']['index']
    output: rma = "{prefix}.rma"
    threads: config['malt']['run']['threads']
    shell:
        "{params.xvfb} malt-run {params.options} -m {params.mode} -i {input.fa} -d {input.index} -o {output.rma}"
