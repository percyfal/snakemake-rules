# -*- snakemake -*-
include: '../ngs.settings'

# Jar program
GATK_JAR_PROGRAM = "GenomeAnalysisTK.jar"
# Subprograms
APPLY_RECALIBRATION = "ApplyRecalibration"
BASE_RECALIBRATOR = "BaseRecalibrator"
CLIP_READS = "ClipReads"
COMBINE_GVCFS = "CombineGVCFs"
COMBINE_VARIANTS = 'CombineVariants'
FASTA_ALTERNATE_REFERENCE_MAKER = 'FastaAlternateReferenceMaker'
GENOTYPE_GVCFS = "GenotypeGVCFs"
HAPLOTYPE_CALLER = 'HaplotypeCaller'
INDEL_REALIGNER = "IndelRealigner"
PRINT_READS = "PrintReads"
READ_BACKED_PHASING = "ReadBackedPhasing"
REALIGNER_TARGET_CREATOR = "RealignerTargetCreator"
SELECT_VARIANTS = "SelectVariants"
UNIFIED_GENOTYPER = "UnifiedGenotyper"
VARIANT_FILTRATION = "VariantFiltration"
VARIANT_EVAL = "VariantEval"
VARIANT_RECALIBRATOR = "VariantRecalibrator"
VARIANTS_TO_TABLE = "VariantsToTable"

config_default = { 
    'gatk' : {
        'bait_regions' : config['ngs.settings']['sequence_capture']['bait_regions'],
        'bam_list' : "",
        'cov_interval' : "regional",
        'dbsnp' : config['ngs.settings']['db']['dbsnp'],
        'home' : os.getenv("GATK_HOME", os.curdir),
        'java_mem' : config['settings']['java']['java_mem'],
        'java_tmpdir' : config['settings']['java']['java_tmpdir'],
        'ref' : config['ngs.settings']['db']['ref'],
        'runtime' : "24:00:00",
        'target_regions' : config['ngs.settings']['sequence_capture']['target_regions'],
        'threads' : config['settings']['threads'],
        'vcfsuffix' : "(.vcf|.vcf.gz)",
    },
}

update_config(config_default, config)
config = config_default

config_default2 = {'gatk': {}}

config_default2 = {
    'gatk' : {
        'jar' : os.path.join(config['gatk']['home'], GATK_JAR_PROGRAM),
        'known_sites' : config['gatk']['dbsnp'],
    },
}

update_config(config_default2, config)
config = config_default2


gatk_args = ["-Xmx" + config['gatk']['java_mem'], "-Djava.io.tmpdir=" + config['gatk']['java_tmpdir'], " "]

if shutil.which("gatk") is None:
    gatk = " ".join(["java"] + gatk_args + ["-jar", join(config['gatk']['jar'])])
else:
    gatk = " ".join(["gatk"] + gatk_args)

config_default3 = {
    'gatk': {
        'cmd' : gatk,
    }
}

update_config(config_default3, config)
config = config_default3

_gatk_config_rule_default = {
    'java_mem' : config['gatk']['java_mem'],
    'options' : '',
    'ref' : config['gatk']['ref'],
    'runtime' : config['gatk']['runtime'],
    'threads' : 1,
}
