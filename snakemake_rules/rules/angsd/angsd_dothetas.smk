# -*- snakemake -*-
include: 'angsd.settings.smk'
include: 'angsd_realsfs.rule'

config_default = {'angsd' :{'dothetas' : _angsd_config_rule_default.copy()}}
config_default['angsd']['dothetas'].update(
    {
        'dothetas' : 1,
        'gl' : config['angsd']['gl'],
	'majorminor' : config['angsd']['majorminor'],
        'dosaf' : config['angsd']['dosaf']['dosaf'],
    })


update_config(config_default, config)
config = config_default

rule angsd_dothetas:
    """angsd dothetas: calculate neutrality tests

    NB: currently assumes input is bam.
    """
    params: cmd = config['angsd']['cmd'],
            options = " ".join(['-doThetas {}'.format(str(config['angsd']['dothetas']['dothetas'])),
                                '-doSaf {}'.format(str(config['angsd']['dothetas']['dosaf'])),
                                '-GL {}'.format(str(config['angsd']['dothetas']['gl'])),
                                config['angsd']['dothetas']['options']]),
            runtime = config['angsd']['dothetas']['runtime']
    input: fofn = "{prefix}.fofn", sfs = "{prefix}.sfs", anc = config['angsd']['anc']
    output: thetas = "{prefix}.thetas.gz"
    threads: config['angsd']['dothetas']['threads']
    conda: "env.yaml"
    shell: "{params.cmd} -bam {input.fofn} -pest {input.sfs} {params.options} -nThreads {threads} -anc {input.anc} -out {wildcards.prefix}"
