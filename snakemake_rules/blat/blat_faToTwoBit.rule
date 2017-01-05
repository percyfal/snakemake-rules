# -*- snakemake -*-
include: "blat.settings"

config_default = {'blat' :{'faToTwoBit' : _blat_config_rule_default.copy()}}
config_default['blat']['faToTwoBit'].update({'cmd' : BLAT_FATOTWOBIT})

update_config(config_default, config)
config = config_default


rule blat_faToTwoBit:
    """Run blat faToTwoBit.

    Convert fasta format sequence files to the dense, randomly
    accessible .2bit format that gfClient can use

    """
    params: cmd = config['blat']['faToTwoBit']['cmd'],
            options = config['blat']['faToTwoBit']['options'],
            runtime = config['blat']['faToTwoBit']['runtime'],
    input: fa = "{prefix}.fa"
    output: twobit = "{prefix}.2bit"
    conda: "env.yaml"
    threads: config['blat']['faToTwoBit']['threads']
    shell: "{params.cmd} {params.options} {input.fa} {output.twobit}"
