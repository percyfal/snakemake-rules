# -*- snakemake -*-
import re
from csv import DictReader

include: "ercc.settings"

config_default = {'ercc' :{'create_ref' : _ercc_config_rule_default.copy()}}
config_default['ercc']['create_ref'].update(
    {
        'ref' : 'ERCC_spikes.fa',
        'retmax' : 100,
    })

update_config(config_default, config)
config = config_default


rule ercc_create_ref:
    """Create ERCC reference file"""
    params: runtime = config['ercc']['create_ref']['runtime']
    input: txt = os.path.basename(config['ercc']['source'])
    output: gb = protected(re.sub("(.fa$|.fasta$)", ".gb", config['ercc']['create_ref']['ref'])),
            fa = protected(config['ercc']['create_ref']['ref'])
    threads: config['ercc']['create_ref']['threads']
    run:
        with open(input.txt, 'r') as fh:
            reader = DictReader(fh, delimiter='\t')
            gb2ercc = {row['GenBank']:row['ERCC_ID'] for row in reader}
            if config['settings']['email'] is None:
                s = "ercc: No email set; must tell NCBI mail address! Either set USER_EMAIL environment variable or configuration variable settings.email"
                logger.warning(s)
                raise
            else:
                from Bio import Entrez, SeqIO
                Entrez.email = config['settings']['email']
                logger.info("ercc: connecting to Entrez; querying for genbank accessions")
                handle = Entrez.efetch(db="nuccore", id=",".join(list(gb2ercc.keys())), rettype="gb", retmode="text")
                records = SeqIO.parse(handle, "gb")
                outrecords = []
                for rec in records:
                    rec.name = gb2ercc[rec.name]
                    rec.id = rec.name
                    outrecords.append(rec)
                SeqIO.write(outrecords, output.gb, "gb")
                SeqIO.write(outrecords, output.fa, "fasta")
