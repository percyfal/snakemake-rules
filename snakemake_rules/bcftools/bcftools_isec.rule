# -*- snakemake -*-
include: "bcftools.settings"

config_default = {'bcftools' :{'isec' : _bcftools_config_rule_default.copy()}}
config_default['bcftools']['isec']['options'] = "-O z"

update_config(config_default, config)
config = config_default

rule bcftools_isec:
    """bcftools isec: Create intersections, unions and complements of VCF files

    """
    params: cmd = config['bcftools']['cmd'],
            options = config['bcftools']['isec']['options'],
            runtime = config['bcftools']['isec']['runtime']
    wildcard_constraints: suffix = "(vcf|vcf.gz|fofn)"
    input: vcf = "{prefix}.{suffix}"
    output: isec = "{prefix}.{suffix}.isec"
    threads: config['bcftools']['isec']['threads']
    conda: "env.yaml"
    shell:
        "if [[ {wildcards.suffix} == \"fofn\" ]]; then\n"
        "command=\"{params.cmd} isec {params.options} $(cat {input.vcf} | tr \"\\n\", \" \") -p {output.isec}\"; eval ${{command}};\n"
        "else\n"
        "{params.cmd} isec {params.options} {input.vcf} -p {output.isec}\n"
        "fi"

