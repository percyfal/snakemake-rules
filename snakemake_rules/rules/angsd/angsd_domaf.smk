# -*- snakemake -*-
include: 'angsd.settings.smk'

config_default = {'angsd' :{'domaf' : _angsd_config_rule_default.copy()}}
config_default['angsd']['domaf'].update(
    {
        'domaf' : 1,
        'gl' : config['angsd']['gl'],
	'majorminor' : config['angsd']['majorminor'],
    })

update_config(config_default, config)
config = config_default

rule angsd_domaf:
    """angsd domaf: estimate the Minor Allele Frequencies.

    A LRT is performed on the MAFs; any maf significantly different
    from zero is indicative of a polymorpic site.

    NB: currently assumes input is bam.
    """
    params: cmd = config['angsd']['cmd'],
            options = " ".join(['-doMaf {}'.format(str(config['angsd']['domaf']['domaf'])),
	    	      	    	'-doMajorMinor {}'.format(str(config['angsd']['domaf']['majorminor'])),
                                '-GL {}'.format(str(config['angsd']['domaf']['gl'])),
                                config['angsd']['domaf']['options']]),
            runtime = config['angsd']['domaf']['runtime']
    input: fofn = "{prefix}.fofn"
    output: pos = "{prefix}.mafs.gz"
    threads: config['angsd']['domaf']['threads']
    conda: "env.yaml"
    shell: "{params.cmd} -bam {input.fofn} {params.options} -nThreads {threads} -out {wildcards.prefix}"
