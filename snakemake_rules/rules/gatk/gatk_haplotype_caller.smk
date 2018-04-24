# -*- snakemake -*-
include: 'gatk.settings.smk'
include: '../picard/picard_build_bam_index.smk'

config_default = {'gatk': {'haplotype_caller' : _gatk_config_rule_default.copy()}}

update_config(config_default, config)
config = config_default

cmd = re.sub("-Xmx[0-9a-zA-Z]+", "-Xmx{mem}".format(mem=config['gatk']['haplotype_caller']['java_mem']), config['gatk']['cmd'])

rule gatk_haplotype_caller:
    """Run GATK HaplotypeCaller"""
    wildcard_constraints:
        mode = "(.g|)",
        suffix = "(.vcf.gz|.vcf)"
    params: cmd = cmd + " -T " + HAPLOTYPE_CALLER,
            options = config['gatk']['haplotype_caller']['options'],
            runtime = config['gatk']['haplotype_caller']['runtime'],
    input: bam = "{prefix}.bam", bai = "{prefix}.bai",
           ref = config['gatk']['haplotype_caller']['ref'],
           refidx = config['gatk']['haplotype_caller']['ref'] + ".fai",
           dict = os.path.splitext(config['gatk']['haplotype_caller']['ref'])[0] + ".dict"
    output: gvcf = "{prefix}.haplotype_caller{mode}{suffix}",
            tbi = "{prefix}.haplotype_caller{mode}{suffix}.tbi"
    threads: config['gatk']['haplotype_caller']['threads']
    conda: "env.yaml"
    log: "{prefix}{mode}.vcf.log"
    shell: "mode=''; if [[ '{wildcards.mode}' == '.g' ]]; then mode='-ERC GVCF'; fi; {params.cmd} {params.options} ${{mode}} -I {input.bam} -R {input.ref} -o {output.gvcf} -log {log}"

