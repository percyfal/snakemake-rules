# -*- snakemake -*-
include: "vsearch.settings.smk"

config_default = {
    'vsearch' : {
        'fastq_filter' : {
            'options' : "",
            'runtime' : "01:00:00",
        },
    },
}

update_config(config_default, config)
config = config_default

rule vsearch_fastq_filter:
    """vsearch --fastq_filter: filter FASTQ file; output is FASTQ file"""
    params: cmd = config['vsearch']['cmd'],
            options = config['vsearch']['fastq_filter']['options'],
            runtime = config['vsearch']['fastq_filter']['runtime']
    threads: config['vsearch']['threads']
    input: fastq = "{prefix}" + config["ngs.settings"]["fastq_suffix"]
    output: fastq = "{prefix}.filter" + re.sub(config["ngs.settings"]["zip_re"], "", config["ngs.settings"]["fastq_suffix"]),
            log = "{prefix}.fastq_filter.txt"
    conda: "env.yaml"
    shell: "{params.cmd} --threads {threads} --fastq_filter {input.fastq} {params.options} --log {output.log} --fastqout {output.fastq}"
