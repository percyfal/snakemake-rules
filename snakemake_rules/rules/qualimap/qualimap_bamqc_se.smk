# -*- snakemake -*-
include: "qualimap.settings.smk"

config_default = {
    'qualimap' : {
        'bamqc' : {
            'options' : [],
            'threads' : 4,
            'runtime' : '01:00:00',
            'outputs' : ['coverage_across_reference', 'coverage_histogram',
                         'duplication_rate_histogram', 'genome_fraction_coverage',
                         'homopolymer_indels', 'mapped_reads_clipping_profile',
                         'mapped_reads_gc-content_distribution', 'mapped_reads_nucleotide_content',
                         'mapping_quality_across_reference', 'mapping_quality_histogram'],
        },
    },
}

update_config(config_default, config)
config = config_default

# Check for -gff flag and -os flag
_input = {'bam': "{prefix}.bam"}
_output = {'outdir': "{prefix}.bam.qualimap",
           'genome_results': "{prefix}.bam.qualimap/genome_results.txt",
           'html': "{prefix}.bam.qualimap/qualimapReport.html"}

m = re.search("(?P<gff>-gff|--feature-file)\s+(?P<featurefile>.*)", " ".join(config['qualimap']['bamqc']['options']))
if m:
    logger.info("{} given in option string; adding feature file output")
    _input.update({'gff': m.group("featurefile")})
if re.search("(?P<os>-os|--outside-stats)", " ".join(config['qualimap']['bamqc']['options'])):
    _output.update({'outside_results': "{prefix}.bam.qualimap/outside_results.txt",
                    'outside_windows': "{prefix}.bam.qualimap/outside_windows.txt",
                    'outside_html': "{prefix}.bam.qualimap/qualimapReportOutsideRegions.html"})

output_map = {
    'coverage_across_reference': "{prefix}.bam.qualimap/raw_data_qualimapReport/coverage_across_reference.txt",
    'coverage_histogram': "{prefix}.bam.qualimap/raw_data_qualimapReport/coverage_histogram.txt",
    'duplication_rate_histogram': "{prefix}.bam.qualimap/raw_data_qualimapReport/duplication_rate_histogram.txt",
    'genome_fraction_coverage': "{prefix}.bam.qualimap/raw_data_qualimapReport/genome_fraction_coverage.txt",
    'homopolymer_indels': "{prefix}.bam.qualimap/raw_data_qualimapReport/homopolymer_indels.txt",
    'mapped_reads_clipping_profile': "{prefix}.bam.qualimap/raw_data_qualimapReport/mapped_reads_clipping_profile.txt",
    'mapped_reads_gc-content_distribution': "{prefix}.bam.qualimap/raw_data_qualimapReport/mapped_reads_gc-content_distribution.txt",
    'mapped_reads_nucleotide_content': "{prefix}.bam.qualimap/raw_data_qualimapReport/mapped_reads_nucleotide_content.txt",
    'mapping_quality_across_reference': "{prefix}.bam.qualimap/raw_data_qualimapReport/mapping_quality_across_reference.txt",
    'mapping_quality_histogram': "{prefix}.bam.qualimap/raw_data_qualimapReport/mapping_quality_histogram.txt",
}

_output.update({k:output_map[k] for k in config['qualimap']['bamqc']['outputs']})

rule qualimap_bamqc_se:
    """Qualimap: run bamqc on bam file"""
    wildcard_constraints:
        end = "(se|pe)"
    params: cmd = config['qualimap']['cmd'],
            options = " ".join(config['qualimap']['bamqc']['options'] +
                               [" --java-mem-size={}".format(available_mem(config['qualimap']['bamqc']['threads'],
                                                                           config['qualimap']['java_mem']))]),
            runtime = config['qualimap']['bamqc']['runtime']
    input: **_input
    output: **_output
    threads: config['qualimap']['bamqc']['threads']
    conda: "env.yaml"
    shell: "unset DISPLAY; {params.cmd} bamqc -bam {input.bam} -nt {threads} {params.options} -outdir {output.outdir}"

