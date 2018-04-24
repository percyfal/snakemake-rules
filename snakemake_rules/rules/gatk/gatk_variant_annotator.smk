# -*- snakemake -*-
include: 'gatk.settings.smk'

config_default = {'gatk' : {'variant_annotator' :  _gatk_config_rule_default.copy()}}
config_default['gatk']['variant_annotator'].update(
    {'annotations': ["QualByDepth", "FisherStrand", "StrandOddsRatio", "RMSMappingQuality", "MappingQualityRankSumTest", "ReadPosRankSumTest", "InbreedingCoeff"],
    'bamfiles': []}
)


update_config(config_default, config)
config = config_default

cmd = re.sub("-Xmx[0-9a-zA-Z]+", "-Xmx{mem}".format(mem=config['gatk']['variant_annotator']['java_mem']), config['gatk']['cmd'])

rule gatk_variant_annotator:
    """Run GATK VariantAnnotator.

    The bam files on which raw variant calling can be set via
    config['gatk']['variant_annotator']['bamfiles'].

    """
    wildcard_constraints:
        suffix = "(.vcf.gz|.vcf)"
    params: cmd = cmd + " -T " + VARIANT_ANNOTATOR,
            options = config['gatk']['variant_annotator']['options'],
            annotations = " ".join("-A {}".format(x) for x in config['gatk']['variant_annotator']['annotations']),
            runtime = config['gatk']['variant_annotator']['runtime']
    input: vcf = "{prefix}{suffix}", ref = config['gatk']['variant_annotator']['ref'],
           bamfiles = config['gatk']['variant_annotator']['bamfiles']
    output: vcf = "{prefix}.annotated{suffix}"
    threads: config['gatk']['variant_annotator']['threads']
    conda: "env.yaml"
    log: "{prefix}.annotated.log"
    shell: "command=\"{params.cmd} {params.options} -nt {threads} -R {input.ref} -V {input.vcf} -o {output} -log {log} {params.annotations} $(echo {input.bamfiles} | sed -e 's/[^ ][^ ]*/-I &/g')\"; eval \"${{command}}\""
