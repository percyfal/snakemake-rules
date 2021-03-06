# -*- snakemake -*-
include: "ucsc.settings.smk"

rule ucsc_link:
    """Link reference to ucsc directory"""
    input: ref = config['ucsc']['ref']
    output: reflink = temp(config['ucsc']['index'])
    run:
        dname = os.path.dirname(output.reflink)
        if not os.path.exists(dname):
            try:
                os.makedirs(dname)
            except OSError:
                if not os.path.isdir(dname):
                    raise
        tgt = os.path.join(os.path.dirname(output.reflink), os.path.basename(input.ref))
        if not os.path.exists(tgt):
            os.symlink(input.ref, tgt)


localrules: ucsc_link            
