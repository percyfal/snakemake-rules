# -*- coding: utf-8 -*-
# Helper functions for extracting information from rules
import os
from os.path import abspath, dirname, join
import re
import sys
import yaml
from snakemake.parser import parse


regex_ft = re.compile("[ ]*(?P<prefix>\{[a-zA-Z_0-9]+\})(?P<ext>[_\/\.a-zA-Z0-9 ]+)[ ]*")

def determine_filetype(io):
    """Determine filetype from pattern.
    
    Params:
      io (str): string representation of input/output

    Returns:
      str: file type extension
    """
    m = regex_ft.search(io)
    # Regular extensions
    if m:
        return m.groupdict().get('ext').lstrip(".")
    else:
        # Could be just prefix string
        return None
        
# Map file extension to file extension fixture; keys are concatenated
# with | and compiled to regular expression
re_filetypes = {
    "bam" : "bam",
    "bai" : "bai",
    "bed" : "bed",
    "fa" : "fasta",
    "fasta" : "fasta",
    "fai" : "fastaindex",
    "fofn" : "fofn",
    "sam" : "sam",
    "samtools_stats.txt" : "samtools_stats",
    "vcf" : "vcf",
    "vcf.gz" : "vcf",
    "tbi" : "tabix",
}


regex = re.compile(r"\.({})\"".format("|".join("{}".format(x) for x in sorted(re_filetypes.keys()))))
regex_name = re.compile(r"[ ]*name=\'(?P<name>[a-zA-Z_0-9]+)[ ]*'")
regex_target_list = re.compile(r"[\"\'][ ]*(?P<target>[^=]+)[ ]*[\"\']")
regex_wildcard_constraints = re.compile(r"(?P<key>[a-zA-z0-9_]+)[ ]*=[ ]*[\"\'](?P<constraint>[^=]+)[ ]*[\"\']")
# The tacit assumption is that @workflow() always ends with a newline
regex_workflow = re.compile(r"@workflow\.(?P<block>[a-zA-z0-9]+)\((?P<content>[\s\S]*?)\)\n")

def parse_rule(rule, prefix=None):
    """Generate information for rule stored in a dictionary.

    Params:
      rule (str): file containing a snakemake rule 
    
    Results:
      input (list): list of input files
      output (list): list of output targets
    """
    d = {}
    codemap = {'rule': "", 'input': "", 'output': "", 'wildcard_constraints': ""}
    rn = os.path.basename(rule).replace(".rule", "")
    app = os.path.basename(os.path.dirname(rule))
    code, linemake, rulecount = parse(rule)

    l = regex_workflow.findall(code)
    for k, v in l:
        codemap[k] = re.sub("[\n\t ]", "", v)

    m_name = regex_name.search(codemap['rule'])
    d['name'] = m_name.group("name")
    d['output'] = regex_target_list.findall(codemap["output"])
    d['input'] = regex_target_list.findall(codemap["input"])
    d['wildcard_constraints'] = {k:v for k, v in regex_wildcard_constraints.findall(codemap["wildcard_constraints"])}

    # Output missing; return input case if possible
    return d


