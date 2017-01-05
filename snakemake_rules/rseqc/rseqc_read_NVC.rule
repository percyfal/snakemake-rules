# -*- snakemake -*-
include: "rseqc.settings"

config_default = {'rseqc' :{'read_NVC' : _rseqc_config_rule_default.copy()}}

update_config(config_default, config)
config = config_default


rule rseqc_read_NVC:
    """Run RSeQC read_NVC.py"""
    params: cmd = READ_NVC,
            options = config['rseqc']['read_NVC']['options'],
            runtime = config['rseqc']['read_NVC']['runtime']
    input: "{prefix}.bam"
    output: xls = "{prefix}_rseqc/read_NVC.NVC.xls"
    threads: config['rseqc']['read_NVC']['threads']
    conda: "env.yaml"
    shell: " {params.cmd} {params.options} -i {input} -o $(dirname {output.xls})/read_NVC"

