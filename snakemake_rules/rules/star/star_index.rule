# -*- snakemake -*-
include: "star.settings"

config_default = {
    'star' : {
        'star_index' : {
            'threads' : config['settings']['threads'],
            'options' : "--genomeSAindexNbases 14",
            'runtime' : config['star']['runtime'],
            'sjdbGTFfile' : config['ngs.settings']['annotation']['transcript_annot_gtf'],
            'sjdbOverhang' : 99,
            'SAname' : "SA",
        },
    },
}

update_config(config_default, config)
config = config_default

rule star_index:
    """Generate STAR genome index. By default will generate index in a
    directory '../rnaseq/star' relative to the directory of the reference
    sequence.

    Remember: for small genomes the parameter --genomeSAindexNbases
    must be adjusted; it is calculated as min(14, log2(GenomeLength)/2 - 1)

    """
    params: cmd = config['star']['cmd'],
            options = " ".join([\
                                config['star']['star_index']['options'],\
                                "--sjdbGTFfile {}".format(config['star']['star_index']['sjdbGTFfile']) if config['star']['star_index']['sjdbGTFfile'] else "",\
                                "--sjdbOverhang {}".format(config['star']['star_index']['sjdbOverhang']) if config['star']['star_index']['sjdbGTFfile'] else ""\
                                ]),
            genomedir = dirname(join(os.curdir, config['star']['index'])),
            runtime = config['star']['star_index']['runtime']
    input: ref = [config['star']['ref']] + ([config['star']['extra_ref']] if config['star']['extra_ref'] else []),
           gtf = config['star']['star_index']['sjdbGTFfile'] if config['star']['star_index']['sjdbGTFfile'] else []
    output: Genome = config['star']['index'], SA=os.path.join(os.path.dirname(config['star']['index']), "SA")
    log: config['star']['index'] + ".log"
    threads: config['star']['star_index']['threads']
    conda: "env.yaml"
    shell: "{params.cmd} --runThreadN {threads} --runMode genomeGenerate --genomeDir {params.genomedir} --genomeFastaFiles {input.ref} {params.options} > {log}"
