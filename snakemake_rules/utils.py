# Copyright (C) 2017 by Per Unneberg
import os
import csv
import ast
import psutil
import subprocess as sp

def get_samples(config, logger):
    if config['settings'].get('sampleinfo') is None:
        logger.info("[snakemake_rules.settings]: no sampleinfo provided; not setting config['_sampleinfo']")
        return
    try:
        ## Get the sample information
        _include = config.get('samples', [])
        _exclude = config.get('ignore_samples', [])
        _ignored = []
        if len(_include) == 0:
            config['samples'] = []
        with open(config['settings']['sampleinfo']) as csvfile:
            dialect = csv.Sniffer().sniff(csvfile.read(1024))
            csvfile.seek(0)
            reader = csv.DictReader(csvfile, dialect=dialect)
            rows = [row for row in reader]
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
        if len(_include) == 0:
            config['samples'] = sorted(list(set(config['samples'])))
            logger.info("[snakemake_rules.settings]: no samples listed in config; added {} samples to config['samples']".format(len(config['samples'])))
        logger.info("[snakemake_rules.settings]: successfully loaded {} entries from '{}'".format(len(config['_sampleinfo']), config['settings']['sampleinfo']))
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
