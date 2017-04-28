# Copyright (C) 2016 by Per Unneberg
# Helper functions for filetypes
import os
import re
import yaml
from os.path import join as pjoin
from os.path import abspath, normpath, dirname

try:
    from pytest_ngsfixtures import ROOT_DIR
except ImportError as e:
    print("\n\n   pytest-ngsfixtures not installed; install with 'conda install -c percyfal pytest-ngsfixtures'\n\n")
    ROOT_DIR = os.curdir

with open(os.path.join(os.path.dirname(__file__), "inputconfig.yaml")) as fh:
    INPUTCONFIG = yaml.load(fh)

# Map file extension to fixture file; keys are concatenated with | and
# compiled to regular expression
# 1. pytest_ngsfixtures
fixture_ext = {
    ".bam" : "data/applications/pe/PUR.HG00731.tiny.sort.bam",
    ".bai" : "data/applications/pe/PUR.HG00731.tiny.sort.bai",
    ".bam.bai" : "data/applications/pe/PUR.HG00731.tiny.sort.bam.bai",
    ".bed" : "data/ref/scaffolds-targets.bed",
    ".bed12" : "data/ref/scaffolds-transcripts.tiny.bed12",
    "chrom.sizes": "data/ref/scaffolds.chrom.sizes",
    ".dict": "data/ref/scaffolds.dict",
    ".fa" : "data/ref/scaffolds.fa",
    #".fasta" : "fasta",
    ".fai" : "data/ref/scaffolds.fa.fai",
    ".fofn" : ["data/applications/pe/tiny.sort.bam.fofn",
               "data/applications/pe/PUR.HG00731.tiny.sort.bam",
               "data/applications/pe/PUR.HG00731.tiny.sort.bai",
               "data/applications/pe/PUR.HG00733.tiny.sort.bam",
               "data/applications/pe/PUR.HG00733.tiny.sort.bai"],
    ".genePred": "data/ref/scaffolds-transcripts.tiny.genePred",
    ".gtf": "data/ref/scaffolds-transcripts.tiny.gtf",
    ".refFlat": "data/ref/scaffolds-transcripts.tiny.refFlat",
    ".interval_list": "data/ref/scaffolds-targets.interval_list",
    # "sam" : "sam",
    # Possibly this should go under rule-specific inputs
    ".stats.txt" : "data/applications/samtools/1.3.1/pe/medium.stats.txt",
    ".vcf" : "data/ref/known.scaffolds.vcf",
    ".vcf.gz" : "data/applications/gatk/3.7/pe/medium.haplotype_caller.vcf.gz",
    ".vcf.gz.tbi" : "data/applications/gatk/3.7/pe/medium.haplotype_caller.vcf.gz.tbi",
    #"tbi" : "tabix",
}
for k, v in fixture_ext.items():
    if isinstance(v, list):
        fixture_ext[k] = [pjoin(ROOT_DIR, x) for x in v]
    else:
        fixture_ext[k] = pjoin(ROOT_DIR, v)

local_ext = {
    'config': normpath(pjoin(dirname(__file__), '../config/config.yaml')),
    'sampleinfo': normpath(pjoin(dirname(__file__), '../config/sampleinfo.csv')),
}

fixture_ext.update(local_ext)

regex = re.compile(r"(?P<ext>({}))['\"]?$".format("|".join("{}".format(x) for x in sorted(fixture_ext.keys()))))


def determine_fixture(io):
    """Determine fixture from pattern.

    Params:
      io (str): string representation of input/output

    Returns:
      str: file fixture
    """
    m = regex.search(io)
    # Regular extensions
    if m:
        return fixture_ext[m.groupdict().get('ext')]
    else:
        # Could be just prefix string
        return None


def set_inputmap(rule):
    """Set inputmap which pairs a wildcard with an output file name

    Params:
      rule(dict): rule info
    """
    def _flatten(d):
        l = d['_list']
        for k, v in d.items():
            if k == "_list":
                continue
            else:
                l.append(v)
        return l
    plist = _flatten(rule['input'])
    inputmap = []
    # Add config if any has config
    if any([re.search("config", x) for x in plist]):
        inputmap = [("config", fixture_ext['config'])]
    if not INPUTCONFIG.get(rule['app'], {}).get(rule['name'], {}).get('ft', None) is None:
        inputmap = inputmap + [(x, fixture_ext[x]) for x in INPUTCONFIG[rule['app']][rule['name']]['ft']]
    return inputmap + [(x, determine_fixture(x)) for x in plist]


def set_output(rule, wildcards):
    """Set outputs.

    Params:
      rule(dict): rule info
      wildcards(dict): dictionary of wildcards
    """
    def _flatten(d):
        l = d['_list']
        for k, v in d.items():
            if k == "_list":
                continue
            else:
                l.append(v)
        return l
    plist = _flatten(rule['output'])
    if not INPUTCONFIG.get(rule['app'], {}).get(rule['name'], {}).get('output', None) is None:
        return INPUTCONFIG[rule['app']][rule['name']]['output']
    return [o.format(**wildcards) for o in plist]
