# -*- snakemake -*-
include: "rseqc.settings.smk"

config_default = {'rseqc' :{'clipping_profile' : _rseqc_config_rule_default.copy()}}

update_config(config_default, config)
config = config_default


rule rseqc_clipping_profile:
    """Run RSeQC clipping_profile.py"""
    params: cmd = CLIPPING_PROFILE,
            options = config['rseqc']['clipping_profile']['options'],
            runtime = config['rseqc']['clipping_profile']['runtime'],
    wildcard_constraints:
        mode = "(SE|PE)"
    input: bam = "{prefix}.bam"
    output: xls = "{prefix}_rseqc/clippingprofile_{mode}.clipping_profile.xls",
    threads: config['rseqc']['clipping_profile']['threads']
    conda: "env.yaml"
    shell: " {params.cmd} {params.options} -i {input.bam} -o $(dirname {output.xls})/clippingprofile_{wildcards.mode} -s {wildcards.mode}"

