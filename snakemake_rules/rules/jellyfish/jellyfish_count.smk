# -*- snakemake -*-
include: "jellyfish.settings.smk"

config_default = {
    'jellyfish': {
        'count': {
            'options': "-m 21 -s 100M -C",
            'runtime' : config['jellyfish']['runtime'],
        },
    },
}

update_config(config_default, config)
config = config_default

re_ext = "(.fq|.fq.gz|.fa|.fa.gz|.fasta|.fastq|.fasta.gz|.fastq.gz|.txt)"
re_ext_fastq = "(.fa$|.fasta$|.fq$|.fastq$)"

rule jellyfish_count:
    params: cmd = config['jellyfish']['cmd'] + " count",
            options = config['jellyfish']['count']['options'],
            runtime = config['jellyfish']['count']['runtime'],
    input: seq = "{prefix}{suffix}"
    output: jf = "{prefix}{suffix," + re_ext + "}.mer_counts.jf"
    threads: config['jellyfish']['threads']
    conda: "env.yaml"
    shell:
        "if file --mime-type -b -L {input.seq} | grep -q gzip; then zcat {input.seq} | {params.cmd} {params.options} -t {threads} /dev/fd/0 -o {output.jf}; "
        "else "
        "if [[ {input.seq} =~ " + re_ext_fastq + " ]]; then {params.cmd} {params.options} -t {threads} {input.seq} -o {output.jf}; "
        "else zcat $( cat {input.seq} | xargs ) | {params.cmd} {params.options} -t {threads} /dev/fd/0 -o {output.jf}; fi; fi"
