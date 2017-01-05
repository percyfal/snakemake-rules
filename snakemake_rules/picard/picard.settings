# -*- snakemake -*-
include: '../ngs.settings'
include: '../comp/comp.settings'

# Jar program names
BED_TO_INTERVAL_LIST = "BedToIntervalList"
BUILD_BAM_INDEX = "BuildBamIndex"
INTERVAL_LIST_TO_BED = "IntervalListToBed"
SORT_SAM = "SortSam"
MERGE_SAM_FILES = "MergeSamFiles"
REORDER_SAM = "ReorderSam"
MARK_DUPLICATES = "MarkDuplicates"
CREATE_SEQUENCE_DICTIONARY = "CreateSequenceDictionary"
COLLECT_INSERT_SIZE_METRICS = "CollectInsertSizeMetrics"
COLLECT_ALIGNMENT_SUMMARY_METRICS = "CollectAlignmentSummaryMetrics"
CALCULATE_HS_METRICS = "CalculateHsMetrics"
ADD_OR_REPLACE_READ_GROUPS = "AddOrReplaceReadGroups"

config_default = { 
    'picard' : {
        'options' : "VALIDATION_STRINGENCY=SILENT",
        'home' : os.getenv("PICARD_HOME", os.curdir),
        'java_mem' : config['settings']['java']['java_mem'],
        'java_tmpdir' : config['settings']['java']['java_tmpdir'],
        'ref' : config['ngs.settings']['db']['ref'],
        'runtime' : "01:00:00",
    },
}

update_config(config_default, config)
config = config_default

picard_args = ["-Xmx" + config['picard']['java_mem'], "-Djava.io.tmpdir=" + config['picard']['java_tmpdir'], " "]

if shutil.which("picard") is None:
    picard = " ".join(["java"] + picard_args + ["-jar", join(config['picard']['home'], "picard.jar ")])
else:
    picard = " ".join(["picard"] + picard_args)
    
config_default2 = {
    'picard' : {
        'cmd' : picard,
    },
}

update_config(config_default2, config)
config = config_default2


_picard_config_rule_default = {
    'options' : '',
    'runtime' : config['picard']['runtime'],
    'threads' : 1,
}
