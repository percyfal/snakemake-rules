# -*- snakemake -*-
include: "cloudbiolinux.settings.smk"

config_default = {
    'cloudbiolinux' : { 
        'annotation' : {
            'transcript_annot_gtf' : None,
        },
    },
}

update_config(config_default, config)
config = config_default

# FIXME: rule in general is wrong
rule cloudbiolinux_update_annotation_gtf:
    input: gtf = [],
           ref = [config['ngs.settings']['db']['ref']] + [x for x in config['ngs.settings']['db']['extra_ref']]
    output: gtf = protected(config['cloudbiolinux']['annotation']['transcript_annot_gtf'])
    run:
        from Bio import SeqIO
        with open(output.gtf, "w") as out:
            with open (input.gtf, "r") as fh:
                for line in fh.readlines():
                    out.write(line)
            for f in input.ref:
                # Use genbank if present - it should be
                gb = re.sub("(.fa$|.fasta$)", ".gb", f)
                if os.path.exists(gb):
                    seqs = SeqIO.parse(open(gb), format = "gb")
                else:
                    logger.warning("cloudbiolinux: No genbank file {} for {}; skipping inclusion in {}".format(gb, f, output.gtf))
                for seq in seqs:
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
                    out.write('{name}\tcloudbiolinux.cloudbiolinux_update_annotation_gtf\t{feature}\t{start}\t{end}\t{score}\t{strand}\t{frame}\tgene_id "{gene_id}"; transcript_id "{transcript_id}";\n'.format(**d))
                    d['feature'] = 'exon'
                    out.write('{name}\tcloudbiolinux.cloudbiolinux_update_annotation_gtf\t{feature}\t{start}\t{end}\t{score}\t{strand}\t{frame}\tgene_id "{gene_id}"; transcript_id "{transcript_id}";\n'.format(**d))
