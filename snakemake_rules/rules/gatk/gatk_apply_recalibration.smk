# -*- snakemake -*-
include: 'gatk.settings.smk'

config_default = {'gatk' :{'apply_recalibration': _gatk_config_rule_default.copy()}}
config_default['gatk']['apply_recalibration'].update(
    {
        'options': "",
    })

update_config(config_default, config)
config = config_default

rule gatk_apply_recalibration:
    """Run GATK ApplyRecalibration."""
    wildcard_constraints:
        mode = "(snp|indel)",
        suffix = "(.vcf|.vcf.gz)"
    params: cmd = config['gatk']['cmd'] + " -T " + APPLY_RECALIBRATION,
            options =  config['gatk']['apply_recalibration']['options'],
            runtime = config['gatk']['apply_recalibration']['runtime']
    input: vcf = "{prefix}.{mode}{suffix}", tbi = "{prefix}.{mode}{suffix}.tbi",
           ref = config['gatk']['apply_recalibration']['ref'],
           fai = config['gatk']['apply_recalibration']['ref'] + ".fai",
           recal = "{prefix}.{mode}.recal", tranches="{prefix}.{mode}.tranches"
    output: vcf="{prefix}.{mode}.recal{suffix}"
    threads: config['gatk']['apply_recalibration']['threads']
    conda: "env.yaml"
    shell: "{params.cmd} -R {input.ref} -input {input.vcf} -recalFile {input.recal} -mode {wildcards.mode} -tranchesFile {input.tranches} {params.options} -o {output.vcf}"

