# -*- snakemake -*-
include: 'gatk.settings'

config_default = {'gatk': {'realigner_target_creator': _gatk_config_rule_default.copy()}}

update_config(config_default, config)
config = config_default


rule gatk_realigner_target_creator:
    """Run GATK RealignerTargetCreator"""
    params: cmd = config['gatk']['cmd'] + " -T " + REALIGNER_TARGET_CREATOR,
            options = " ".join(["-R", config['gatk']['realigner_target_creator']['ref'],
            config['gatk']['realigner_target_creator']['options']]),
            runtime = config['gatk']['realigner_target_creator']['runtime']
    input: bam="{prefix}.bam",
           bai="{prefix}.bai",
           d=config['gatk']['realigner_target_creator']['ref'].replace(".fa", ".dict"),
           fai="{}.fai".format(config['gatk']['realigner_target_creator']['ref'])
    output: intervals = "{prefix}.intervals"
    threads: config['gatk']['realigner_target_creator']['threads']
    conda: "env.yaml"
    shell: "{params.cmd} {params.options} -I {input.bam} -o {output.intervals}"
