# -*- snakemake -*-
include: "star.settings"
include: "star_index.rule"

config_default = { 
    'star' : {
        'align' : {
            'threads': config['settings']['threads'],
            'options': "",
            'outSAMtype': "Unsorted",
            'runtime' : config['star']['runtime'],
        },
    },
}

update_config(config_default, config)
config = config_default

def _star_suffix():
    if config['star']['align']['outSAMtype'].find("Unsorted") > -1:
        return ".Aligned.out.bam"
    elif config['star']['align']['outSAMtype'].find("SortedByCoordinate") > -1:
        return ".Aligned.sortedByCoord.out.bam"

rule star_align_pe:
    """Run STAR alignment"""
    params: cmd = config['star']['cmd'],
            genomedir = dirname(join(os.curdir, config['star']['index'])),
            options = " ".join([config['star']['align']['options'],
                                "--readFilesCommand {cmd}".format(cmd=config['comp']['compression']['prog_map'].get(os.path.splitext(config['ngs.settings']['fastq_suffix'])[1], "")),
                                "--outSAMtype BAM {}".format(config['star']['align']['outSAMtype']),
                                "--quantMode TranscriptomeSAM" if config['ngs.settings']['rnaseq']['_transcriptome_quantification'] else "",
                                ]),
            runtime = config['star']['align']['runtime']
    input: read1 = "{prefix}" + config['ngs.settings']['read1_label'] + config['ngs.settings']['fastq_suffix'],
           read2 = "{prefix}" + config['ngs.settings']['read2_label'] + config['ngs.settings']['fastq_suffix'],
           index = config['star']['index']
    output: bam = ["{prefix}" + _star_suffix()],
            txbam = ["{prefix}" + ".Aligned.toTranscriptome.out.bam"] if config['ngs.settings']['rnaseq']['_transcriptome_quantification'] else [],
            log = "{prefix}.Log.final.out"
    log: "{prefix}.Log.final.out"
    threads: config['star']['align']['threads']
    conda: "env.yaml"
    shell: "{params.cmd} --runThreadN {threads} --genomeDir {params.genomedir} --readFilesIn {input.read1} {input.read2} {params.options} --outFileNamePrefix {wildcards.prefix}."

