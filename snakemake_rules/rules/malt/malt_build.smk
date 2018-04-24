# -*- snakemake -*-
include: 'malt.settings'

config_default = {'malt' :{'build' : _malt_config_rule_default.copy()}}
config_default['malt']['build']['sequenceType'] = "Protein"

update_config(config_default, config)
config = config_default

rule malt_build:
    """Build an index for MALT (MEGAN alignment tool)

    The input is a list of fastq files (possibly zipped) that need to
    be set in the configuration (config['malt']['ref']).

    Output is a malt index, also defined in the configuration file
    (config['malt']['index']).

    """
    params: runtime = config['malt']['build']['runtime'],
            xvfb = "" if config['malt']['X11'] else "xvfb-run -a ",
            options = config['malt']['build']['options'],
            sequenceType = config['malt']['build']['sequenceType'],
            indexdir = config['malt']['index']
    input: fa = config['malt']['ref']
    output: ref = expand(os.path.join(config['malt']['index'], "ref.{ext}"), ext=["db", "idx", "inf"]),
            taxonomy = expand(os.path.join(config['malt']['index'], "taxonomy.{ext}"), ext=["map", "idx", "tre"]),
            index = expand(os.path.join(config['malt']['index'], "index{num}.idx"), num=[0, 1, 2, 3]),
            table = expand(os.path.join(config['malt']['index'], "table{num}.{ext}"), num=[0, 1, 2, 3], ext=["db", "idx"]),
    threads: config['malt']['build']['threads']
    shell:
        "{params.xvfb} malt-build {params.options} -t {threads} -s {params.sequenceType} -i {input.fa} -d {params.indexdir}"
