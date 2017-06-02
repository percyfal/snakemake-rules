# -*- coding: utf-8 -*-
import os
import re
import csv
import ast
import subprocess as sp
import psutil
import math
from snakemake import available_cpu_count

BYTE = 1024


def get_samples(config, logger):
    if config['settings'].get('sampleinfo') is None:
        logger.info("[snakemake_rules.settings]: no sampleinfo provided; not setting config['_sampleinfo']")
        return
    try:
        # Get the sample information
        _include = config.get('samples', [])
        if _include is None:
            _include = []
        assert isinstance(_include, list), \
            print("config['samples'] should be a list")
        _exclude = config.get('ignore_samples', [])
        if _exclude is None:
            _exclude = []
        assert isinstance(_exclude, list), \
            print("config['ignore_samples] should be a list")
        _ignored = []
        if len(_include) == 0:
            config['samples'] = []
        with open(config['settings']['sampleinfo']) as csvfile:
            dialect = csv.Sniffer().sniff(csvfile.read(10240),
                                          delimiters=[',', ';'])
            # save for later use by QC rules
            config['_sampleinfo_delim'] = dialect.delimiter
            csvfile.seek(0)
            reader = csv.DictReader(csvfile, dialect=dialect)
            rows = [row for row in reader]
            # check if sample in config but not in sampleinfo:
            if len(_include)>0 and len([s for s in _include if s not in [row['SM'] for row in rows]])>0:
                logger.info("[snakemake_rules.settings]: warning: sample listed in config['samples'] not in '{}'".format(config['settings']['sampleinfo']))
            if len(_exclude)>0 and len([s for s in _exclude if s not in [row['SM'] for row in rows]])>0:
                logger.info("[snakemake_rules.settings]: warning: sample listed in config['ignore_samples'] not in '{}'".format(config['settings']['sampleinfo']))
            # check overlap include exclude:
            if len(_include)>0 and len(_exclude)>0 and not set(_include).isdisjoint(_exclude):
                logger.info("[snakemake_rules.settings]: warning: same sample(s) listed in config['samples'] and config['ignore_samples']")
            config['_sampleinfo'] = []
            for row in rows:
                # Only excluded samples listed
                if len(_include) == 0 and row['SM'] not in _exclude:
                    config['_sampleinfo'].append(row)
                    config['samples'].append(row['SM'])
                # include takes precedence over exclude
                elif row['SM'] in _include:
                    config['_sampleinfo'].append(row)
                elif row['SM'] in _exclude:
                    _ignored.append(row['SM'])
                else:
                    logger.warning("[snakemake_rules.settings]: {} neither in include or exclude section of sampleinfo file".format(row['SM']))
            for s in _ignored:
                logger.info("[snakemake_rules.settings]: ignoring sample {} in current run".format(s))
        logger.info("[snakemake_rules.settings]: successfully loaded {} sample(s) from '{}'".format(len(config['_sampleinfo']), config['settings']['sampleinfo']))
        if len(_include) == 0:
            config['samples'] = sorted(list(set(config['samples'])))
            logger.info("[snakemake_rules.settings]: no samples listed in config; added {} samples to config['samples']".format(len(config['samples'])))
        if len(_ignored) > 0:
            logger.info("[snakemake_rules.settings]: ignored {} sample(s) listed in '{}'".format(len(_ignored), config['settings']['sampleinfo']))        
    except Exception as e:
        logger.info("[snakemake_rules.settings]: parsing sampleinfo failed; not setting config['_sampleinfo']: {}".format(e))


def python2_path(config, logger):
    """Add python 2 path if possible"""
    try:
        output = sp.check_output(['conda', 'env', 'list', '--json'])
        envs = ast.literal_eval(output.decode("utf-8"))
        py2 = [x for x  in envs['envs'] if re.search("{}{}$".format(os.sep, config['settings']['conda']['python2']), x)][0]
        py2bin = os.path.join(py2, "bin")
        os.environ["PATH"] = ":".join([os.environ["PATH"], py2bin])
        logger.info("[snakemake_rules.settings]: Add python2 bin path '{}' to PATH".format(py2bin))
    except:
        logger.info("[snakemake_rules.settings]: failed to add conda python2 environment '{py2}' to PATH".format(py2=config['settings']['conda']['python2']))


# From https://github.com/giampaolo/psutil/blob/master/scripts/meminfo.py
def bytes2human(n):
    # http://code.activestate.com/recipes/578019
    # >>> bytes2human(10000)
    # '9.8K'
    # >>> bytes2human(100001221)
    # '95.4M'
    symbols = ('K', 'M', 'G', 'T', 'P', 'E', 'Z', 'Y')
    prefix = {}
    for i, s in enumerate(symbols):
        prefix[s] = 1 << (i + 1) * 10
    for s in reversed(symbols):
        if n >= prefix[s]:
            value = float(n) / prefix[s]
            return '%.1f%s' % (value, s)
    return "%sB" % n


def mem_per_core(prefix="G"):
    """Calculate available memory per core"""
    conversion = 1e6
    if prefix == "G" or prefix == "g":
        conversion = 1e6
    elif prefix == "M" or prefix == "m":
        conversion = 1e3
    else:
        pass
    cores = psutil.cpu_count()
    mem = psutil.virtual_memory().total
    return math.floor(mem / BYTE / conversion / cores)


def available_mem(cores, mem, fmtstring=True):
    """Calculate available memory for a process

    Params:
      cores (int): number of cores
      mem (str): set memory as string with conversion (M, G, g)
    fmtstring (bool): return memory as formatted string
    """
    prefix = "G"
    m = re.match("[0-9]+([a-zA-Z]*)", str(mem))
    if m:
        prefix = m.groups()[0]
    requested_mem_per_core = int(re.sub("[a-zA-Z]*", "", str(mem)))
    core_mem = mem_per_core(prefix)
    requested_cores = min(cores, available_cpu_count())
    mem = min(requested_cores * core_mem,
              requested_cores * requested_mem_per_core)
    if fmtstring:
        return "{}{}".format(mem, prefix)
    else:
        return mem
