# -*- snakemake -*-
include: "bcftools.settings"

config_default = {'bcftools' :{'isec_AB' : _bcftools_config_rule_default.copy()}}

update_config(config_default, config)
config = config_default


rule bcftools_isec_AB:
    """bcftools isec_AB: Create intersections, unions and complements of VCF files

    """
    wildcard_constraints:
        sampleA = "(" + "|".join(config["samples"]) + ")",
        sampleB = "(" + "|".join(config["samples"]) + ")"
    params: cmd = config['bcftools']['cmd'],
            options = config['bcftools']['isec_AB']['options'],
            runtime = config['bcftools']['isec_AB']['runtime']
    input: A = "{prefix}{sampleA}.vcf.gz", Atbi = "{prefix}{sampleA}.vcf.gz.tbi",
           B = "{prefix}{sampleB}.vcf.gz", Btbi = "{prefix}{sampleB}.vcf.gz.tbi",
    output: isec_0000 = "{prefix}{sampleA}_{sampleB}.isec/0000.vcf.gz",
            isec_0001 = "{prefix}{sampleA}_{sampleB}.isec/0001.vcf.gz",
            isec_0002 = "{prefix}{sampleA}_{sampleB}.isec/0002.vcf.gz",
            isec_0003 = "{prefix}{sampleA}_{sampleB}.isec/0003.vcf.gz",
    threads: config['bcftools']['isec_AB']['threads']
    conda: "env.yaml"
    shell:
        "{params.cmd} isec -O z {params.options} {input.A} {input.B} -p $(dirname {output.isec_0000})"
