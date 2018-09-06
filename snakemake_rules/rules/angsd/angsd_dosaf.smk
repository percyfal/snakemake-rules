# -*- snakemake -*-
include: 'angsd.settings.smk'

config_default = {'angsd' :{'dosaf' : _angsd_config_rule_default.copy()}}
config_default['angsd']['dosaf'].update(
    {
        'dosaf' : 1,
        'gl' : config['angsd']['gl'],
    })

update_config(config_default, config)
config = config_default

rule angsd_dosaf:
    """angsd dosaf: estimate the SFS.

    NB: currently assumes input is bam.
    """
    params: cmd = config['angsd']['cmd'],
            options = " ".join(['-dosaf', str(config['angsd']['dosaf']['dosaf']),
                                '-GL', str(config['angsd']['dosaf']['gl']),
                                config['angsd']['dosaf']['options']]),
            runtime = config['angsd']['dosaf']['runtime']
    threads: config['angsd']['dosaf']['threads']
    input: fofn = "{prefix}.fofn", anc = config['angsd']['anc']
    output: pos = "{prefix}.saf.pos.gz", saf = "{prefix}.saf.gz", idx = "{prefix}.saf.idx"
    conda: "env.yaml"
    shell: "{params.cmd} -bam {input.fofn} -anc {input.anc} {params.options} -out {wildcards.prefix}"
