# -*- snakemake -*-
include: 'gatk.settings'
include: '../picard/picard_create_sequence_dictionary.rule'

config_default = {'gatk': {'fasta_alternate_reference_maker' : _gatk_config_rule_default.copy()}}

update_config(config_default, config)
config = config_default


rule gatk_fasta_alternate_reference_maker:
    """Run GATK FastaAlternateReferenceMaker"""
    wildcard_constraints: vcfsuffix = "(.vcf|.vcf.gz)"
    params: cmd = config['gatk']['cmd'] + " -T " + FASTA_ALTERNATE_REFERENCE_MAKER,
            options = " ".join([config['gatk']['fasta_alternate_reference_maker']['options']]),
            runtime = config['gatk']['fasta_alternate_reference_maker']['runtime']
    input: vcf = "{prefix}{vcfsuffix}",
           dict = os.path.splitext(config['gatk']['fasta_alternate_reference_maker']['ref'])[0] + ".dict",
           ref = config['gatk']['fasta_alternate_reference_maker']['ref']
    output: "{prefix}.fa"
    threads: config['gatk']['fasta_alternate_reference_maker']['threads']
    conda: "env.yaml"
    shell: "{params.cmd} {params.options} -R {input.ref} -V {input.vcf} -o {output}"

