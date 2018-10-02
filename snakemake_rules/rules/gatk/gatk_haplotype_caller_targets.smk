# -*- snakemake -*-
include: "gatk.settings.smk"
include: "../pybedtools/pybedtools_make_bed_targets.smk"
include: "../samtools/samtools_faidx.smk"
include: '../picard/picard_build_bam_index.smk'

config_default = {'gatk': {'haplotype_caller_targets' : _gatk_config_rule_default.copy()}}

update_config(config_default, config)
config = config_default

cmd = re.sub("-Xmx[0-9a-zA-Z]+", "-Xmx{mem}".format(mem=config['gatk']['haplotype_caller_targets']['java_mem']), config['gatk']['cmd'])

rule gatk_haplotype_caller_targets:
    """Run GATK HaplotypeCaller on a bam file given targets in a bed file

    Takes as input a single bam file, a reference, the reference
    index, and a bed file containing target regions in which variant
    calling is performed. The targets wildcard is often a number,
    resulting from the partitioning of a larger bed files into several
    bed files of subregions.

    
    """
    wildcard_constraints:
        mode = "(.g|)",
        targets = "\d+",
        suffix = "(.vcf.gz|.vcf)"
    params: cmd = cmd + " -T " + HAPLOTYPE_CALLER,
            options = config['gatk']['haplotype_caller_targets']['options'],
            runtime = config['gatk']['haplotype_caller_targets']['runtime'],
    input: bam = "{prefix}.bam",
           bai = "{prefix}.bai",
           ref = config['gatk']['haplotype_caller_targets']['ref'],
           fai = config['gatk']['haplotype_caller_targets']['ref'] + ".fai",
           targets = "{prefix}.{targets}.bed"
    output: vcf = "{prefix}.haplotype_caller.{targets}{mode}{suffix}",
            tbi = "{prefix}.haplotype_caller.{targets}{mode}{suffix}.tbi"
    threads: config['gatk']['haplotype_caller_targets']['threads']
    conda: "env.yaml"
    log: "{prefix}.haplotype_caller.{targets}{mode}.vcf.log"
    shell: "mode=''; if [[ '{wildcards.mode}' == '.g' ]]; then mode='-ERC GVCF'; fi; {params.cmd} {params.options} ${{mode}} -I {input.bam} -L {input.targets} -R {input.ref} -o {output.vcf} -log {log}"
