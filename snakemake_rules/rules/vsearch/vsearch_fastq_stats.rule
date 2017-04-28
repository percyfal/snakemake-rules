# -*- snakemake -*-
include: "vsearch.settings"

config_default = {
    'vsearch' : {
        'fastq_stats' : {
            'options' : "",
            'runtime' : "01:00:00",
        },
    },
}

update_config(config_default, config)
config = config_default

rule vsearch_fastq_stats:
    """vsearch --fastq_stats: FASTQ quality statistics"""
    params: cmd = config['vsearch']['cmd'],
            options = config['vsearch']['fastq_stats']['options'],
            runtime = config['vsearch']['fastq_stats']['runtime']
    threads: config['vsearch']['threads']
    input: fastq = "{prefix}" + config["ngs.settings"]["fastq_suffix"]
    output: log = "{prefix}.fastq_stats.txt"
    conda: "env.yaml"
    shell: "{params.cmd} --threads {threads} --fastq_stats {input.fastq} {params.options} --log {output.log}"
