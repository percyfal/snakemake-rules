# -*- snakemake -*-
include: "rseqc.settings"

config_default = {'rseqc' :{'read_GC' : _rseqc_config_rule_default.copy()}}

update_config(config_default, config)
config = config_default


rule rseqc_read_GC:
    """Run RSeQC read_GC.py"""
    params: cmd = READ_GC,
            options = config['rseqc']['read_GC']['options'],
            runtime = config['rseqc']['read_GC']['runtime'],
    input: "{prefix}.bam"
    output: xls = "{prefix}_rseqc/read_GC.GC.xls"
    threads: config['rseqc']['read_GC']['threads']
    conda: "env.yaml"
    shell: " {params.cmd} {params.options} -i {input} -o $(dirname {output.xls})/read_GC"

