# -*- snakemake -*-
# See http://sandberg.cmb.ki.se/multo/
#
include: '../main'
include: '../comp/comp.settings.smk'

config_default = {
    'bio.ngs.tools.multo' : {
        'home' : "",
        'cmd'  : "MULTo1.0.py",
        'options' : " -t refGene -s Mmusculus",
        'threads' : config['settings']['threads'],
    },
}

update_config(config_default, config)
config = config_default

rule multo_transcript_level:
    params: cmd = os.path.join(config['multo']['home'], 'src', config['multo']['cmd']),
            options = config['multo']['options']
    threads: config['multo']['threads']
    output: outdir = os.path.join("{outdir}", "{assembly}_{kmin}-{kmax}")
    shell: "{params.cmd} {params.options} -T -a {wildcards.assembly} -k {wildcards.kmin} -m {wildcards.kmax} -p {threads} -o {output.outdir}"
