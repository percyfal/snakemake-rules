# -*- snakemake -*-
include: "tuxedo.settings.smk"

config_default = {
    'tuxedo' : {
        'tophat' : {
            'options' : '',
            'threads' : config['settings']['threads'],
            'output_dir' : os.path.join(os.curdir, "tophat_out"),
        },
    },
}

update_config(config_default, config)
config = config_default


rule tuxedo_tophat:
    """Run tophat paired end alignment.

    NB: requires python2.
    """
    params: cmd = 'tophat2' if config['tuxedo']['version2'] else 'tophat',
            options = config['tuxedo']['tophat']['options'],
            index = str(config['tuxedo']['index']),
            outdir = config['tuxedo']['tophat']['output_dir']
    input: read1 = "{prefix}" + config['ngs.settings']['read1_label'] + config['ngs.settings']['fastq_suffix'],\
           read2 = "{prefix}" + config['ngs.settings']['read2_label'] + config['ngs.settings']['fastq_suffix'],\
           index = str(config['tuxedo']['index']) + ".1.bt2" if config['tuxedo']['version2'] else str(config['tuxedo']['index']) + ".1.ebwt"
    output: tophatdir = "{prefix}.tophat2" if config['tuxedo']['version2'] else "{prefix}.tophat",
            hits = os.path.join("{prefix}.tophat2" if config['tuxedo']['version2'] else "{prefix}.tophat", "accepted_hits.bam")
    run: 
        rg = config['tuxedo']['rg_fn'](os.path.basename(wildcards.prefix))
        shell("{cmd} {opt} {rg} -o {out}.tmp {index} {read1} {read2} 2> {out}.log && mv {out}.tmp/* {out} && mv {out}.log {out}".format(cmd=params.cmd, opt=params.options, rg=rg, out=output.tophatdir, prefix=wildcards.prefix, index=params.index, read1=input.read1, read2=input.read2))

