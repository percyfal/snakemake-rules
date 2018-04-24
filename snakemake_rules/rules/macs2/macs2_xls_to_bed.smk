# -*- snakemake -*-
include: "macs2.settings.smk"
include: "macs2_callpeak_treatment_only_bam_pe.rule"

config_default = {'macs2' :{'xls_to_bed' : _macs2_config_rule_default.copy()}}

update_config(config_default, config)
config = config_default


rule macs2_xls_to_bed:
    """MACS: convert peak xls file to bed.

    NB: this rule has output set to temporary by default
    """
    params: runtime = config['macs2']['xls_to_bed']['runtime']
    input: xls = "{prefix}.xls"
    output: bed = temporary("{prefix}.bed")
    threads: config['macs2']['xls_to_bed']['threads']
    conda: "env.yaml"
    shell: 'cat {input.xls} | grep -v "#" | grep -v "abs_summit" | awk "BEGIN {{OFS=\\"\\t\\"}}; /^\s*$/ {{next; }};  {{print \$1, \$2, \$3, \$10, \\"*\\", \\"*\\"}}" > {output.bed}'
        
    
