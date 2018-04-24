# -*- snakemake -*-
# See https://groups.google.com/forum/#!topic/rna-star/tvajn49WTYk for
# setting up RSEM with STAR alignments
# Take into account different aligners
include: "rsem.settings.smk"

aligner = config['ngs.settings'].get('aligner', 'star')
align_section = '' + aligner
align_include = os.path.join(os.pardir, aligner + ".rules")

include: align_include

rule rsem_calculate_expression_bowtie:
    """Run RSEM on bowtie output. Requires a fifo hack to work: see
    http://atgcio.blogspot.se/2013/08/fifos-and-mapping-with-bowtie-using.html
    """
    params: cmd = config['rsem']['calculate-expression']['cmd'],
            options = config['rsem']['calculate-expression']['options'],
            index = config['rsem']['index'] + config['rsem']['ref_sfx'],
            runtime = config['rsem']['runtime']
    input: read1 = "{prefix}" + config['ngs.settings']['read1_label'] + config['ngs.settings']['fastq_suffix'],\ 
           read2 = "{prefix}" + config['ngs.settings']['read2_label'] + config['ngs.settings']['fastq_suffix'],\
           index = config['rsem']['index'] + config['rsem']['ref_sfx'],
           bowtie_index = expand(config['rsem']['index'] + "{ext}", ext=config[align_section].get('build', {}).get('ext_v1', ""))
    output: "{prefix}.rsem"
    threads: config['rsem']["threads"]
    conda: "env.yaml"
    run: 
      fifo1 = input.read1 + ".read1.fifo"
      fifo2 = input.read2 + ".read2.fifo"
      shell("rm -f " + fifo1)
      shell("rm -f " + fifo2)
      shell("mkfifo " + fifo1)
      shell("mkfifo " + fifo2)
      shell("zcat " + input.read1 + " > " + fifo1 + " &")
      shell("zcat " + input.read2 + " > " + fifo2 + " &")
      shell("{cmd} {options} -p {threads} --paired-end {fifo1} {fifo2} {index} {sample} &> {log}.tmp && mv {log}.tmp {log}.rsem".format(cmd=params.cmd, options=params.options, threads=threads, fifo1=fifo1, fifo2=fifo2, index=config['rsem']['index'], sample=wildcards.prefix, log=wildcards.prefix))
      shell("rm -f " + fifo1)
      shell("rm -f " + fifo2)
