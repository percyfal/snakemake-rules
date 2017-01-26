# Copyright (C) 2016 by Per Unneberg
import os
from os.path import abspath, dirname, join, isdir, basename
import sys
import re
import logging
import pytest
import shutil
import subprocess as sp
import ast
import yaml
from snakemake.parser import parse

# Add helper module
sys.path.append(os.path.join(os.path.dirname(__file__), 'helpers'))

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

TESTDIR = abspath(dirname(__file__))
RULEDIR = join(abspath(dirname(__file__)), os.pardir, "snakemake_rules")

# Add test source path to pythonpath
sys.path.insert(0, join(TESTDIR, os.pardir))

# Autogenerate all application directories
blacklist = ['__pycache__', 'bio', 'comp']
applications = {x:join(RULEDIR, x) for x in os.listdir(RULEDIR) if isdir(join(RULEDIR, x)) and not x in blacklist}
rules = {k: [join(RULEDIR, k, x) for x in os.listdir(v) if x.endswith(".rule")] for k,v in applications.items() }
with open(os.path.join(TESTDIR, "rules2target.yaml")) as fh:
    rules2targets = yaml.load(fh)

def make_output(rule, prefix="s1"):
    """Generate output file for rule. Assume target is based on prefix
    's1'; else, use rules2targets dictionary"""

    rn = os.path.basename(rule).replace(".rule", "")
    app = os.path.basename(os.path.dirname(rule))
    target = rules2targets.get(app, {}).get(rn, None)
    if target:
        return target
    code, linemake, rulecount = parse(rule)
    m = re.search("@workflow.output\(\s+(?P<output>.*)", code)
    if m is None:
        # return input case
        m = re.search("@workflow.input\(\s+(?P<input>.*)", code)
        output = m.group("input")
    else:
        output = m.group("output")
    m = re.search("\"[ ]*(?P<prefix>\{[a-zA-Z_0-9]+\})+(?P<ext>[_\/\.a-zA-Z0-9 ]+)\"", output)
    # Regular extension; use first one
    if m:
        return "{prefix}{ext}".format(prefix=prefix, ext=m.group("ext"))
    # expand case; skip for now
    m = re.search("expand", output)
    if m:
        return None
    # Config case
    m = re.search("[a-zA-Z =]*(?P<config>config[^\)]+)", output)
    if m:
        return "config"
    return None


# Add option to run slow tests; by default these are turned off
def pytest_addoption(parser):
    parser.addoption("--slow", action="store_true",
                     help="run slow tests", dest="slow")
    parser.addoption("--slow-only", action="store_true",
                     help="run slow tests only", dest="slow_only")
    parser.addoption("-V", "--show-workflow-output", action="store_true",
                     help="show workflow output",
                     dest="show_workflow_output")
    parser.addoption("-P", "--python2-conda", action="store",
                     default="py2.7",
                     help="name of python2 conda environment [default: py2.7]",
                     dest="python2_conda")
    parser.addoption("-A", "--application", action="store",
                     default=False,
                     help="application to test",
                     dest="application")
    parser.addoption("-T", "--threads", action="store",
                     default="1",
                     help="number of threads to use",
                     dest="threads")
    parser.addoption("-R", "--rule", action="store",
                     default=False,
                     help="run a specific rule",
                     dest="rule")


def pytest_runtest_setup(item):
    """Skip tests if they are marked as slow and --slow is not given"""
    if getattr(item.obj, 'slow', None) and not (item.config.getvalue('slow') or item.config.getvalue('slow_only')):
        pytest.skip('slow tests not requested')
    if not getattr(item.obj, 'slow', None) and item.config.getvalue('slow_only'):
        pytest.skip('run only slow tests')


def pytest_report_header(config):
    try:
        output = sp.check_output(['conda', 'env', 'list', '--json'])
        envs = ast.literal_eval(output.decode("utf-8"))
        py2 = [x for x  in envs['envs'] if re.search("{}{}$".format(os.sep, config.getoption("--python2-conda")), x)][0]
        py2bin = join(py2, "bin")
        os.environ["PATH"] = ":".join([os.environ["PATH"], py2bin])
        py2str = "python2: {}".format(py2bin)

    except:
        py2str = "WARNING: No conda python2 environment found! Workflow tests depending on python2 programs will fail!"
    return "\n".join([py2str])

# Add namespace with test files and director
def pytest_namespace():
    return {
        'testdir': TESTDIR,
        'ruledir': RULEDIR,
        'snakefile': join(TESTDIR, "Snakefile"),
        'snakefile_regions' : join(TESTDIR, "Snakefile_regions"),
        'config' : join(TESTDIR, "config.yaml"),
        'config_regions' : join(TESTDIR, "config_regions.yaml"),
        'rules' : rules,
        'make_output' : make_output,
    }

##################################################
# Setup test fixtures
##################################################
d = {}
for f in os.listdir(abspath(join(dirname(__file__), "data"))):
    d[f] = abspath(join(dirname(__file__), "data", f))

copyfiles = ['s1.revealjs.Rmd']
skipfiles = ['Snakefile']
    
d.update({
    "s1.bam" : d['s1.rg.sort.bam'],
    "s1.bai" : d['s1.rg.sort.bai'],
    "s1.fofn" : d['bamfiles.fofn'],
    "s1.bam.fofn" : d['bamfiles.fofn'],
    "s1.bdg" : d['s1.rg.sort.bedGraph'],
    "s1.bedGraph" : d['s1.rg.sort.bedGraph'],
    "s1.bed" : d['ref.bed'],
    "s1.dict" : d['ref.dict'],
    "s1.gtf" : d['ref-transcripts.gtf'],
    "s1.genePred" : d['ref-transcripts.genePred'],
    "s1.fasta" : d['ref.fa'],
    "s1.fa" : d['ref.fa'],
    "s1.fastq.gz" : d['s1_1.fastq.gz'],
    "s1.fofn" : d['bamfiles.fofn'],
    "s1.g.vcf.fofn" : d['s1.vcf.fofn'],
    "s1.interval_list" : d['ref.interval_list'],
    "s1.sam" : d['s1.rg.sort.sam'],
    "s1.sort.bam" : d['s1.rg.sort.bam'],
    "s1.sort.bai" : d['s1.rg.sort.bai'],
    "s1.sort.bam.bai" : d['s1.rg.sort.bai'],
    "s1.vcf" : d['s1.rg.sort.vcf'],
    "s1.vcf.gz" : d['s1.rg.sort.g.vcf.gz'],
    "s1.wig" : d['s1.rg.sort.wig'],
})


##############################
# sample
##############################
@pytest.fixture(scope="function", autouse=False)
def data(tmpdir_factory):
    """Setup input data

    FIXME: currently a test directory is setup for *every* atomic
    test. However, setting it up only once will fail parallel tests.

    """
    p = tmpdir_factory.mktemp('data')

    for k, v in d.items():
        if basename(v) in skipfiles:
            continue
        if basename(v) in copyfiles:
            shutil.copyfile(v, str(p.join(k)))
        else:
            p.join(k).mksymlinkto(v)
    return p

##############################
# Snakefile data
##############################
@pytest.fixture(scope="function", autouse=False)
def snakefile_data(tmpdir_factory):
    """Setup input data for snakefiles"""
    p = tmpdir_factory.mktemp('snakefile_data')
    data = abspath(join(dirname(__file__), "data"))
    p.join("ref.fa").mksymlinkto(join(data, "ref.fa"))
    p.join("s1_1.fastq.gz").mksymlinkto(join(data, "s1_1.fastq.gz"))
    p.join("s1_2.fastq.gz").mksymlinkto(join(data, "s1_2.fastq.gz"))
    
    path = abspath(join(dirname(__file__), "examples"))
    config = join(path, "config.yaml")
    snakefile = join(path, "Snakefile")
    p.join("config2.yaml").mksymlinkto(config)
    p.join("Snakefile2").mksymlinkto(snakefile)

    config = join(path, "config_regions.yaml")
    snakefile = join(path, "Snakefile_regions")
    p.join("config1.yaml").mksymlinkto(config)
    p.join("Snakefile1").mksymlinkto(snakefile)

    return p
