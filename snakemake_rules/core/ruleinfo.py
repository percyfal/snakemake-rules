# -*- coding: utf-8 -*-
# Helper functions for extracting information from rules
import os
from os.path import abspath, dirname, join
import re
import sys
import yaml
from snakemake.parser import parse

regex_name = re.compile(r"[ ]*name=\'(?P<name>[a-zA-Z_0-9]+)[ ]*'")
regex_target_list = re.compile(r"(?P<key>[a-zA-Z0-9_]+=)?(?P<target>[^,]+)")
regex_wildcard_constraints = re.compile(r"(?P<key>[a-zA-Z0-9_]+)[ ]*=[ ]*[\"\'](?P<constraint>[^=]+)[ ]*[\"\']")
# The tacit assumption is that @workflow() always ends with a newline
regex_workflow = re.compile(r"@workflow\.(?P<block>[a-zA-z0-9]+)\((?P<content>[\s\S]*?)\)\n")

def get_targets(targets):
    d = {'_list': []}
    l = regex_target_list.findall(targets)
    for k, v in l:
        if k == "":
            d['_list'].append(v)
        else:
            d[k.rstrip("=")] = v
    return d


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
        codemap[k] = re.sub("[\t\n ]", "", v)

    m_name = regex_name.search(codemap['rule'])
    d['name'] = m_name.group("name")
    d['output'] = get_targets(codemap["output"])
    d['input'] = get_targets(codemap["input"])
    d['wildcard_constraints'] = {k:v for k, v in regex_wildcard_constraints.findall(codemap["wildcard_constraints"])}
    if d['wildcard_constraints'] == '':
        d['wildcard_constraints'] = {}

    # Output missing; return input case if possible
    return d


