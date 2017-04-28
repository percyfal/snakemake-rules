# -*- snakemake -*-
include: 'gatk.settings'

config_default = {'gatk' : {'base_recalibrator' : _gatk_config_rule_default.copy()}}
config_default['gatk']['base_recalibrator'].update({
                'options' : " ".join([
                "-L {target}".format(target=config['gatk']['target_regions']) if not config['gatk']['target_regions'] == "" else "",
                "-knownSites {known}".format(known=config['gatk']['known_sites'] if not config['gatk']['known_sites'] == "" else ""),
            ]),
})


update_config(config_default, config)
config = config_default


rule gatk_base_recalibrator:
    """Run GATK BaseRecalibrator"""
    params: cmd = config['gatk']['cmd'] + " -T " + BASE_RECALIBRATOR,
            options = " ".join(["-R", config['gatk']['base_recalibrator']['ref'],
            config['gatk']['base_recalibrator']['options']]),
            runtime = config['gatk']['base_recalibrator']['runtime']
    threads: config['gatk']['base_recalibrator']['threads']
    conda: "env.yaml"
    input: bam = "{prefix}.bam", bai = "{prefix}.bai",
           known = config['gatk']['known_sites'], knowntbi = config['gatk']['known_sites'] + ".tbi"
    output: grp = "{prefix}.recal_data.grp"
    shell: "{params.cmd} {params.options} -I {input.bam} -o {output.grp}"
