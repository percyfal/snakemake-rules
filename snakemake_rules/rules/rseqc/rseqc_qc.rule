# -*- snakemake -*-
include: "rseqc.settings"

RULES = {
    'rseqc_clipping_profile': (True, "{prefix}_rseqc/clippingprofile.clipping_profile.xls"),
    'rseqc_geneBody_coverage': (True, "{prefix}_rseqc/geneBody_coverage.geneBodyCoverage.txt"),
    'rseqc_junction_annotation': (True, "{prefix}_rseqc/junction_annotation_refseq.txt"),
    'rseqc_read_GC': (True, "{prefix}_rseqc/read_GC.GC.xls"),
    'rseqc_read_NVC': (True, "{prefix}_rseqc/read_NVC.NVC.xls"),
    'rseqc_read_distribution': (True, "{prefix}_rseqc/read_distribution.txt"),
    'rseqc_read_duplication': (True, "{prefix}_rseqc/read_dup.pos.DupRate.xls"),
    'rseqc_read_quality': (True, "{prefix}_rseqc/read_quality.qual.r"),
}


config_default = {
    'rseqc' : {
        'rseqc_qc': {k:v[0] for k,v in RULES.items() if v[0]},
        'transcriptome_only': False,
    },
}

update_config(config_default, config)
config = config_default

# Remove rules for which qc make no sense, like geneBody coverage
if config['rseqc']['transcriptome_only']:
    for r in ['rseqc_read_distribution', 'rseqc_geneBody_coverage', 'rseqc_junction_annotation']: # ['rseqc_tin']
        config['rseqc'][r] = False

        
for rule in RULES.keys():
    if not RULES[rule][0]:
        continue
    if not rule in workflow._rules:
        include: rule + ".rule"

        
rule rseqc_qc:
    """Run selection of RSeQC commands on a bam file"""
    input: **{k:RULES[k][1] for k,v in config['rseqc']['rseqc_qc'].items() if v}
    output: "{prefix}_rseqc/rseqc_qc.txt"
    shell: "echo `date` > {output}"

localrules: rseqc_qc
