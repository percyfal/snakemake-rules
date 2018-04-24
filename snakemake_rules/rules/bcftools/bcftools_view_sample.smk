# -*- snakemake -*-
include: "bcftools.settings"

config_default = {'bcftools' :{'view_sample' : _bcftools_config_rule_default.copy()}}

update_config(config_default, config)
config = config_default

rule bcftools_view_sample:
    """bcftools view sample: View a sample in variant file

    The sample of choice is given by the {sample} tag; types can be
    selected via the {type} tag.

    """
    params: cmd = config['bcftools']['cmd'],
            options = config['bcftools']['view_sample']['options'],
            runtime = config['bcftools']['view_sample']['runtime'],
    wildcard_constraints: suffix = "(vcf|vcf.gz)",
    			  type = "(.|.snps.|.indels.|.mnps.|.ref.|.bnd.|.other.)"
    input: vcf = "{prefix}.{suffix}"
    output: vcf = "{prefix}.{sample}{type}{suffix}"
    threads: config['bcftools']['view_sample']['threads']
    conda: "env.yaml"
    shell:
        "if [[ {wildcards.suffix} == \"vcf.gz\" ]]; then\n"
        "outputtype=\"z\";\n"
        "else\n"
        "outputtype=\"v\";\n"
        "fi\n"
        "if [[ {wildcards.type} == \".\" ]]; then\n"
        "{params.cmd} view {params.options} -O $outputtype -s {wildcards.sample} {input.vcf} -o {output.vcf}\n"

        "else\n"
        "{params.cmd} view {params.options} -O $outputtype -v {wildcards.type} -s {wildcards.sample} {input.vcf} -o {output.vcf}\n"
        "fi"

