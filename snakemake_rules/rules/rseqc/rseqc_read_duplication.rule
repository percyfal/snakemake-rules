# -*- snakemake -*-
include: "rseqc.settings"

config_default = {'rseqc' :{'read_duplication' : _rseqc_config_rule_default.copy()}}

update_config(config_default, config)
config = config_default


rule rseqc_read_duplication:
    """Run RSeQC read_duplication.py"""
    params: cmd = READ_DUPLICATION,
            options = config['rseqc']['read_duplication']['options'],
            runtime = config['rseqc']['read_duplication']['runtime'],
    input: bam = "{prefix}.bam"
    output: pos = "{prefix}_rseqc/read_dup.pos.DupRate.xls",
            seq = "{prefix}_rseqc/read_dup.seq.DupRate.xls",
    threads: config['rseqc']['read_duplication']['threads']
    conda: "env.yaml"
    shell: " {params.cmd} {params.options} -i {input.bam} -o $(dirname {output.seq})/read_dup"
