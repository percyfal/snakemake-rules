# -*- snakemake -*-
include: "rseqc.settings"

config_default = {'rseqc' :{'read_distribution' : _rseqc_config_rule_default.copy()}}

update_config(config_default, config)
config = config_default


rule rseqc_read_distribution:
    """Run RSeQC read_distribution.py
    
    NB: Requries reference gene model in bed format. Also memory
    intensive; needs more resources (approx 17GB).
    """
    params: cmd = READ_DISTRIBUTION,
            options = config['rseqc']['read_distribution']['options'],
            runtime = config['rseqc']['read_distribution']['runtime'],
    input: bam = "{prefix}.bam", refgene = config['rseqc']['refgene']
    output: "{prefix}_rseqc/read_distribution.txt"
    threads: config['rseqc']['read_distribution']['threads']
    conda: "env.yaml"
    shell: " {params.cmd} {params.options} -i {input.bam} -r {input.refgene} &> {output}"

