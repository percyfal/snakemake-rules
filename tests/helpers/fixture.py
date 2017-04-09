# Copyright (C) 2016 by Per Unneberg
# Helper functions for filetypes
import os
import re
from os.path import join as pjoin
from os.path import abspath, normpath

try:
    from pytest_ngsfixtures import ROOT_DIR
except ImportError as e:
    print("\n\n   pytest-ngsfixtures not installed; install with 'conda install -c percyfal pytest-ngsfixtures'\n\n")
    ROOT_DIR = os.curdir

# Map file extension to fixture file; keys are concatenated with | and
# compiled to regular expression
# 1. pytest_ngsfixtures
fixture_ext = {
    ".bam" : "data/applications/pe/PUR.HG00731.tiny.sort.bam",
    ".bai" : "data/applications/pe/PUR.HG00731.tiny.sort.bai",
    ".bam.bai" : "data/applications/pe/PUR.HG00731.tiny.sort.bam.bai",
    ".bed" : "data/ref/ref-targets.bed",
    ".bed12" : "data/ref/ref-transcripts.tiny.bed12",
    "chrom.sizes": "data/ref/ref.chrom.sizes",
    ".dict": "data/ref/ref.dict",
    ".fa" : "data/ref/ref.fa",
    #".fasta" : "fasta",
    ".fai" : "data/ref/ref.fa.fai",
    ".fofn" : ["data/applications/pe/tiny.sort.bam.fofn",
               "data/applications/pe/PUR.HG00731.tiny.sort.bam",
               "data/applications/pe/PUR.HG00731.tiny.sort.bai",
               "data/applications/pe/PUR.HG00733.tiny.sort.bam",
               "data/applications/pe/PUR.HG00733.tiny.sort.bai"],
    ".genePred": "data/ref/ref-transcripts.tiny.genePred",
    ".gtf": "data/ref/ref-transcripts.tiny.gtf",
    ".refFlat": "data/ref/ref-transcripts.tiny.refFlat",
    ".interval_list": "data/ref/ref-targets.interval_list",
    # "sam" : "sam",
    # Possibly this should go under rule-specific inputs
    ".stats.txt" : "data/applications/samtools/1.3.1/pe/medium.stats.txt",
    ".vcf" : "data/ref/known-ref.vcf",
    ".vcf.gz" : "data/applications/gatk/3.7/pe/medium.haplotype_caller.vcf.gz",
    #"tbi" : "tabix",
}
for k, v in fixture_ext.items():
    if isinstance(v, list):
        fixture_ext[k] = [pjoin(ROOT_DIR, x) for x in v]
    else:
        fixture_ext[k] = pjoin(ROOT_DIR, v)

local_ext = {
    'config': normpath(abspath('../config/config.yaml')),
}

fixture_ext.update(local_ext)

regex = re.compile(r"(?P<ext>({}))".format("|".join("{}$".format(x) for x in sorted(fixture_ext.keys()))))


def determine_fixture(io):
    """Determine fixture from pattern.

    Params:
      io (str): string representation of input/output

    Returns:
      str: file fixture
    """
    # config case
    if not re.search("config", io) is None:
        return fixture_ext["config"]
    m = regex.search(io)
    # Regular extensions
    if m:
        return fixture_ext[m.groupdict().get('ext')]
    else:
        # Could be just prefix string
        return None
