# -*- snakemake -*-
include: "utils.settings"

rule dbutils_make_transcript_annot_gtf:
    """Database utilities: make transcript annotation file from entries
defined in config['ngs.settings']['annotation']['sources']. Source
files can be in gtf or genbank format.

"""
    input: sources = config['ngs.settings']['annotation']['sources']
    output: gtf = protected(config['ngs.settings']['annotation']['transcript_annot_gtf'])
    run:
        from Bio import SeqIO
        def _gb_to_gtf(seq, out):
            d = {
                'name' : seq.name,
                'score' : 0,
                'strand' : '+',
                'frame' : 0,
                'gene_id' : seq.name,
                'transcript_id' : seq.name,
                'feature' : 'transcript',
                'start' : 1,
                'end' : len(seq),
            }
            if len(seq.features) > 0:
                d['start'] = max(seq.features[0].location.start, 1)
                d['end'] = seq.features[0].location.end
            out.write('{name}\tutils.dbutils_make_transcript_annot_gtf\t{feature}\t{start}\t{end}\t{score}\t{strand}\t{frame}\tgene_id "{gene_id}"; transcript_id "{transcript_id}";\n'.format(**d))
            d['feature'] = 'exon'
            out.write('{name}\tutils.dbutils_make_transcript_annot_gtf\t{feature}\t{start}\t{end}\t{score}\t{strand}\t{frame}\tgene_id "{gene_id}"; transcript_id "{transcript_id}";\n'.format(**d))
            
        with open(output.gtf, "w") as out:
            for f in input.sources:
                # Check for input format
                ext = os.path.splitext(f)[1]
                if ext == ".gtf":
                    logger.info("utils: Treating input source {} as gtf format".format(f))
                    with open (f, "r") as fh:
                        for line in fh.readlines():
                            out.write(line)
                elif ext in ['.gb', '.genbank']:
                    logger.info("utils: Treating input source {} as genbank format".format(f))
                    seqs = SeqIO.parse(open(f), format = "gb")
                    for seq in seqs:
                        _gb_to_gtf(seq, out)
                else:
                    logger.warning("utils: Input source extension {ext} unknown for {f}: skipping inclusion in annotation file".format(ext=ext, f=f))
                    raise
