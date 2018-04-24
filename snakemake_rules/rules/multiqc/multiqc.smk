# -*- snakemake -*-
include: "multiqc.settings.smk"

rule multiqc:
    """multiqc: run multiqc in a directory

    By default runs in the working directory.
    """
    params: cmd = config['multiqc']['cmd'],
            options = config['multiqc']['options'],
            runtime = config['multiqc']['runtime']
    threads: config['multiqc']['threads']
    input: config['multiqc']['inputs']
    output: html = os.path.join("multiqc", "multiqc_report.html")
    conda: "env.yaml"
    shell:
        "{params.cmd} {params.options} . -o $(dirname {output.html})"
