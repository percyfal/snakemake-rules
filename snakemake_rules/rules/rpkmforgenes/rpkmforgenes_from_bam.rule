# -*- snakemake -*-
# rpkmforgenes.py -readcount -fulltranscript -mRNAnorm -rmnameoverlap -bothendsceil -n 15_115 -i unique.bam -p 1 -a refGene_140508_norandom_egfp.txt -u MULTo1.0/mm10_20-255/ -o test_refseq_rpkms.txt
include: "rpkmforgenes.settings"

rule rpkmforgenes_from_bam:
    """Run rpkmforgenes from bam input"""
    params: cmd = config['rpkmforgenes']['cmd'],
            options = " ".join([
                config['rpkmforgenes']['options'],
                " -u {}".format(config['rpkmforgenes']['unique']) if config['rpkmforgenes']['unique'] else "",
                " -a {}".format(config['rpkmforgenes']['annotation']) if config['rpkmforgenes']['annotation'] else "",
            ]),
            runtime = config['rpkmforgenes']['runtime'],
    input: unique = [config['rpkmforgenes']['unique']] if not config['rpkmforgenes']['unique'] is None else [],
           annotation = [config['rpkmforgenes']['annotation']] if config['rpkmforgenes']['annotation'] else [],
           bam = "{prefix}.bam"
    output: rpkmforgenes = "{prefix}.rpkmforgenes"
    log: "{prefix}.rpkmforgenes.log"
    threads: config['rpkmforgenes']['threads']
    conda: "env.yaml"
    shell:
        " ".join(["{params.cmd} {params.options} ",
                  "-bamu -i {input.bam} -o {output.rpkmforgenes} ",
                  " &> {log}"])
