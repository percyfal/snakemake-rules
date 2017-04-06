# -*- coding: utf-8 -*-
# Helper functions for extracting information from rules
import os
from os.path import abspath, dirname, join
import re
import sys
import yaml
from snakemake.parser import parse


regex_ft = re.compile("\"[ ]*(?P<prefix>\{[a-zA-Z_0-9]+\})+(?P<ext>[_\/\.a-zA-Z0-9 ]+)\"[ ]*")

def _determine_filetypes(io):
    """Determine filetype from pattern.
    
    Params:
      io (str): string representation of input/output
    """
    m = regex_ft.findall(io)
    # Regular extensions
    if m:
        return [x[1].lstrip(".") for x in m]
    else:
        return [None]
        
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


regex = re.compile("\.({})\"".format("|".join("{}".format(x) for x in sorted(re_filetypes.keys()))))
regex_input = re.compile("@workflow\.input\(\s+(?P<input>.*)")
regex_output = re.compile("@workflow.output\(\s+(?P<output>.*)")
regex_output_list = re.compile("[\"\'][ ]*(?P<target>\{[a-zA-Z_0-9]+\}[\._\/a-zA-Z0-9]+)[\"\']")

def parse_rule(rule, prefix=None):
    """Generate input/output information for rule.

    Params:
      rule (str): file containing a snakemake rule 
    
    Results:
      input (list): list of input filetypes
      output (list): list of output targets
    """
    rn = os.path.basename(rule).replace(".rule", "")
    app = os.path.basename(os.path.dirname(rule))
    code, linemake, rulecount = parse(rule)

    # Regular input/output
    m_out = regex_output.search(code)
    if not m_out is None:
        output = regex_output_list.findall(m_out.group("output"))
    m_in = regex_input.search(code)
    if not m_in is None: 
        input = _determine_filetypes(m_in.group("input"))

    # Output missing; return input case if possible
    if m_out is None and not m_in is None:
        output = []
    # else:
    #     if m_out:
    #         # Return 
    
    return input, output
        # if m is None:
        # return input case
    #     m = regex_input.search(code)
    #     output = m.group("input")
    # else:
    #     output = m.group("output")
    # m = regex.findall(input)
    # if m:
    #     print(m)
    #     #print(minput.groups())
    # else:
    #     print("no match")
    # m = re.search("\"[ ]*(?P<prefix>\{[a-zA-Z_0-9]+\})+(?P<ext>[_\/\.a-zA-Z0-9 ]+)\"", output)
    # # Regular extension; use first one
    # if m:
    #     return "{prefix}{ext}".format(prefix=prefix, ext=m.group("ext"))
    # # expand case; skip for now
    # m = re.search("expand", output)
    # if m:
    #     return None, None
    # Config case
    # m = re.search("[a-zA-Z =]*(?P<config>config[^\)]+)", output)
    # if m:
    #     return "config"
    return None, None

