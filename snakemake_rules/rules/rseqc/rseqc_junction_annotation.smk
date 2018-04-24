# -*- snakemake -*-
include: "rseqc.settings"

config_default = {'rseqc' :{'junction_annotation' : _rseqc_config_rule_default.copy()}}

update_config(config_default, config)
config = config_default


rule rseqc_junction_annotation:
    """Run RSeQC junction_annotation.py
    NB: Requries reference gene model in bed format
    """
    params: cmd = JUNCTION_ANNOTATION,
            options = config['rseqc']['junction_annotation']['options'],
            runtime = config['rseqc']['junction_annotation']['runtime'],
    input: bam = "{prefix}.bam", refgene = config['rseqc']['refgene']
    output: xls = "{prefix}_rseqc/junction_annotation_refseq.junction.xls",
            txt = "{prefix}_rseqc/junction_annotation_refseq.txt"
    threads: config['rseqc']['junction_annotation']['threads']
    conda: "env.yaml"
    shell: " {params.cmd} {params.options} -i {input.bam} -o $(dirname {output.txt})/junction_annotation_refseq -r {input.refgene} 2> {output.txt}"

