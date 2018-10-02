# -*- snakemake -*-
include: "picard.settings.smk"

config_default = {
    'picard' : {
        'qcrules' : ['picard_collect_insert_size_metrics',
                     'picard_collect_alignment_summary_metrics',
                     'picard_mark_duplicates'],
    },
}

update_config(config_default, config)
config = config_default


for rule in config['picard']['qcrules']:
    if not rule in workflow._rules:
        include: rule + ".smk"

rule picard_do_qc:
    """Run picard metrics commands on a bam file"""
    params: runtime = "01:00:00"
    input: ["{{prefix}}{sfx}".format(sfx=workflow._rules[x].params.suffix) for x in config['picard']['qcrules']]
    output: out="{prefix}.picardqc.txt"
    conda: "env.yaml"
    threads: 1
    run:
        with open (output.out, "w") as fh:
            fh.write("Completed rule 'picard_do_qc' using qcrules {qcrules}".format(qcrules=",".join(config['picard']['qcrules'])))

localrules: picard_do_qc
            
