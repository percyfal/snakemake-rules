# -*- snakemake -*-
include: 'gatk.settings.smk'
include: 'gatk_realigner_target_creator.smk'

config_default = {'gatk' : {'indel_realigner' : _gatk_config_rule_default.copy()}}
config_default['gatk']['indel_realigner'].update({'options' : " ".join(["-L {target}".format(target=config['gatk']['target_regions']) if not config['gatk']['target_regions'] == "" else ""])})

update_config(config_default, config)
config = config_default


rule gatk_indel_realigner:
    """Run GATK IndelRealigner"""
    params: cmd = config['gatk']['cmd'] + " -T " + INDEL_REALIGNER,
            options = " ".join(["-R", config['gatk']['indel_realigner']['ref'],
            config['gatk']['indel_realigner']['options']]),
            runtime  = config['gatk']['indel_realigner']['runtime']
    input: bam = "{prefix}.bam",
           bai = "{prefix}.bai",
           intervals = "{prefix}.intervals"
    output: bam = "{prefix}.realign.bam", bai = "{prefix}.realign.bai"
    threads: config['gatk']['indel_realigner']['threads']
    conda: "env.yaml"
    shell: "{params.cmd} {params.options} -o {output.bam} --targetIntervals {input.intervals} -I {input.bam}"
