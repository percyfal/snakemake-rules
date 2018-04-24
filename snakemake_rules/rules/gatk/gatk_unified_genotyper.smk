# -*- snakemake -*-
include: 'gatk.settings.smk'
include: join(os.pardir, "samtools", "samtools_faidx.smk")
include: join(os.pardir, "picard", "picard_create_sequence_dictionary.smk")
include: join(os.pardir, "picard", "picard_build_bam_index.smk")

config_default = {'gatk' : {'unified_genotyper' : _gatk_config_rule_default.copy()}}
config_default['gatk']['unified_genotyper'].update(
    {
        'options' : " ".join(["-stand_call_conf 30.0 -stand_emit_conf 10.0  --downsample_to_coverage 200 --output_mode EMIT_VARIANTS_ONLY -glm BOTH",
                              "--dbsnp {dbsnp}".format(dbsnp=config['gatk']['dbsnp']) if not config['gatk']['dbsnp'] == "" else "",
                              "-L {target}".format(target=config['gatk']['target_regions']) if not config['gatk']['target_regions'] == "" else "",
            ]),
        })

update_config(config_default, config)
config = config_default

cmd = re.sub("-Xmx[0-9a-zA-Z]+", "-Xmx{mem}".format(mem=config['gatk']['unified_genotyper']['java_mem']), config['gatk']['cmd'])


rule gatk_unified_genotyper:
    """Run GATK UnifiedGenotyper"""
    wildcard_constraints:
        suffix = "(.vcf|.vcf.gz)"
    params: cmd = cmd + " -T " + UNIFIED_GENOTYPER,
            options = config['gatk']['unified_genotyper']['options'],
            runtime = config['gatk']['unified_genotyper']['runtime']
    input: bam = "{prefix}.bam", bai = "{prefix}.bai",
           ref = config['gatk']['unified_genotyper']['ref'],
           fai = config['gatk']['unified_genotyper']['ref'] + '.fai',
           dict = os.path.splitext(config['gatk']['unified_genotyper']['ref'])[0] + '.dict',
    output: vcf = "{prefix}{suffix}"
    threads: config['gatk']['unified_genotyper']['threads']
    conda: "env.yaml"
    shell: "{params.cmd} {params.options} -R {input.ref} -I {input.bam} -o {output.vcf}"

