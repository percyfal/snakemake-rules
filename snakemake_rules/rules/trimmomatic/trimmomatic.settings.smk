# -*- snakemake -*-
include: "../ngs.settings.smk"

TRIMMOMATIC_JAR_PROGRAM = "trimmomatic.jar"

config_default = {
    'bio.ngs.qc.trimmomatic' : {
        'home' : os.getenv("TRIMMOMATIC_HOME", os.curdir),
    	'options' : "-phred33",
        'java_mem' : config['settings']['java']['java_mem'],
        'java_tmpdir' : config['settings']['java']['java_tmpdir'],
        'processing_options' : "LEADING:15 TRAILING:15 MINLEN:36",
    },
}

update_config(config_default, config)
config = config_default

config_default2 = {
    'bio.ngs.qc.trimmomatic' : {
        'jar' : os.path.join(config['bio.ngs.qc.sequenceprocessing']['trimmomatic']['home'], TRIMMOMATIC_JAR_PROGRAM),
    }}

update_config(config_default2, config)
config = config_default

config_default3 = {
    'bio.ngs.qc.trimmomatic' : {
        'cmd' : "java -Xmx" + config['bio.ngs.qc.sequenceprocessing']['trimmomatic']['java_mem'] + " -Djava.io.tmpdir=" + config['bio.ngs.qc.sequenceprocessing']['trimmomatic']['java_tmpdir'] +  " -jar " + config['bio.ngs.qc.sequenceprocessing']['trimmomatic']['jar'],
    }}

update_config(config_default3, config)
config = config_default
