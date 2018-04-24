# -*- snakemake -*-
# Miscallaneous tools for conversion between file formats
from snakemake.utils import R

include: 'utils.settings.smk'

config_default = { 
    'bio.ngs.tools.conversion' : {
    },
}

update_config(config_default, config)
config = config_default

    
rule gtf_to_bed12:
    """Convert gtf file to bed12 format.
    """
    input: gtf = "{prefix}.gtf"
    output: bed12 = "{prefix}.bed12"
    run:
        R("""
        library(GenomicFeatures)
        library(rtracklayer)
        txdb = makeTxDbFromGFF("{input.gtf}", format="gtf")
        tx = asBED(exonsBy(txdb, use.names=TRUE))
        tx.df = do.call("rbind", lapply(tx, function(x){{cbind(as.character(seqnames(x)), start(x), end(x), values(x)$name, 0, as.character(strand(x)), 0, 0, "255,0,0", length(x$blocks), paste(width(x$blocks)[[1]], collapse=","), paste(start(x$blocks)[[1]], collapse=","))}}))
        write.table(tx.df, file="{output.bed12}", quote=FALSE, row.names=FALSE, col.names=FALSE)
        """)

