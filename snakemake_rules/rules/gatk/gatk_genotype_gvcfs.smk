# -*- snakemake -*-
include: 'gatk.settings.smk'
include: '../bcftools/bcftools.settings.smk'

config_default = {'gatk' : {'genotype_gvcfs' :  _gatk_config_rule_default.copy()}}

update_config(config_default, config)
config = config_default

cmd = re.sub("-Xmx[0-9a-zA-Z]+", "-Xmx{mem}".format(mem=config['gatk']['genotype_gvcfs']['java_mem']), config['gatk']['cmd'])

rule gatk_genotype_gvcfs:
    """Run GATK GenotypeGVCFs"""
    wildcard_constraints:
        prefix = ".+[^\.g]$"
    params: cmd = cmd + " -T " + GENOTYPE_GVCFS,
            options = config['gatk']['genotype_gvcfs']['options'],
            runtime = config['gatk']['genotype_gvcfs']['runtime'],
            bcftools = config['bcftools']['cmd']
    input: fofn = "{prefix}.g.vcf.fofn", ref = config['gatk']['genotype_gvcfs']['ref'],
           dict = re.sub("(\.fasta$|\.fa$)", ".dict", config['gatk']['genotype_gvcfs']['ref'])
    output: vcf = "{prefix}.gatk.vcf.gz", tbi = "{prefix}.gatk.vcf.gz.tbi"
    threads: config['gatk']['genotype_gvcfs']['threads']
    conda: "env.yaml"
    log: "{prefix}.vcf.log"
    shell: "{params.cmd} {params.options} -nt {threads} -R {input.ref} -o {output.vcf} -log {log} $(cat {input.fofn} | sed -e \"s/^\(.*\)/--variant \\1/g\" | tr \"\\n\", \" \")"
