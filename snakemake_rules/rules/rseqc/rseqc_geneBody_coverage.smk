# -*- snakemake -*-
include: "rseqc.settings.smk"

config_default = {'rseqc' :{'geneBody_coverage' : _rseqc_config_rule_default.copy()}}

update_config(config_default, config)
config = config_default


rule rseqc_geneBody_coverage:
    """Run RSeQC geneBody_coverage.py. 

    NB: Requries reference gene model in bed format (-r parameter).
    Requires indexed bam, which must have the suffix .bam.bai.


    Memory intensive; needs more resources.

    """
    params: cmd = GENEBODY_COVERAGE,
            options = config['rseqc']['geneBody_coverage']['options'],
            runtime = config['rseqc']['geneBody_coverage']['runtime'],
    input: bam = "{prefix}.bam", bai = "{prefix}.bam.bai", refgene = config['rseqc']['refgene']
    output: txt = "{prefix}_rseqc/geneBody_coverage.geneBodyCoverage.txt"
    threads: config['rseqc']['geneBody_coverage']['threads']
    conda: "env.yaml"
    shell: " {params.cmd} {params.options} -i {input.bam} -o $(dirname {output.txt})/geneBody_coverage -r {input.refgene}"

