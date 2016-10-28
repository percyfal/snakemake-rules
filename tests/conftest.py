# Copyright (C) 2016 by Per Unneberg
import os
from os.path import abspath, dirname, join, isdir
import re
import logging
import pytest
import shutil
import subprocess as sp
import ast

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

TESTDIR = abspath(dirname(__file__))
RULEDIR = join(abspath(dirname(__file__)), os.pardir, "snakemake_rules")
SNAKEFILE = join(TESTDIR, "Snakefile")
SNAKEFILE_REGIONS = join(TESTDIR, "Snakefile_regions")
CONFIG = join(TESTDIR, "config.yaml")
CONFIG_REGIONS = join(TESTDIR, "config_regions.yaml")

# Autogenerate all application directories
blacklist = ['__pycache__', 'bio', 'comp']
applications = {x:join(RULEDIR, x) for x in os.listdir(RULEDIR) if isdir(join(RULEDIR, x)) and not x in blacklist}
rules = {k: [join(RULEDIR, k, x) for x in os.listdir(v) if x.endswith(".rule")] for k,v in applications.items() }


# Add namespace with test files and director
def pytest_namespace():
    return {
        'testdir': TESTDIR,
        'ruledir': RULEDIR,
        'snakefile': SNAKEFILE,
        'snakefile_regions' : SNAKEFILE_REGIONS,
        'config' : CONFIG,
        'config_regions' : CONFIG_REGIONS,
        'rules' : rules,
    }


# Add option to run slow tests; by default these are turned off
def pytest_addoption(parser):
    parser.addoption("--slow", action="store_true",
        help="run slow tests")
    parser.addoption("-W", "--show-workflow-output", action="store_true",
                     help="show workflow output",
                     dest="show_workflow_output")
    parser.addoption("-P", "--python2-conda", action="store",
                     default="py2.7",
                     help="name of python2 conda environment [default: py2.7]",
                     dest="python2_conda")


def pytest_runtest_setup(item):
    """Skip tests if they are marked as slow and --slow is not given"""
    if getattr(item.obj, 'slow', None) and not item.config.getvalue('slow'):
        pytest.skip('slow tests not requested')


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


##################################################
# Setup test fixtures
##################################################
# Input files
chr11 = abspath(join(dirname(__file__), "data", "chr11.fa"))
ref_transcripts = abspath(join(dirname(__file__), "data", "ref-transcripts.gtf"))
sample1_1 = abspath(join(dirname(__file__), "data", "s1_1.fastq.gz"))
sample1_2 = abspath(join(dirname(__file__), "data", "s1_2.fastq.gz"))
sample2_1 = abspath(join(dirname(__file__), "data", "s2_1.fastq.gz"))
sample2_2 = abspath(join(dirname(__file__), "data", "s2_2.fastq.gz"))
sampleinfo = abspath(join(dirname(__file__), "data", "sampleinfo.csv"))

##############################
# Reference data
##############################
@pytest.fixture(scope="session", autouse=True)
def ref(tmpdir_factory):
    """Setup reference data

    Reference data layout:

    ref/
    - chr11.fa
    - ref-transcripts.gtf
    """
    p = tmpdir_factory.mktemp('ref', numbered=False)
    p.join("chr11.fa").mksymlinkto(chr11)
    p.join("ref-transcripts.gtf").mksymlinkto(ref_transcripts)
    if pytest.config.getoption("-v"):
        logger.info("Created reference data folder\n\t" + "\n\t".join(str(x) for x in sorted(p.visit())))
    return p


##############################
# sample
##############################
@pytest.fixture(scope="session", autouse=True)
def data(tmpdir_factory):
    """
    Setup input data
    """
    p = tmpdir_factory.mktemp('data')
    p.join("sampleinfo.csv").mksymlinkto(sampleinfo)
    p.join("s1_1.fastq.gz").mksymlinkto(sample1_1)
    p.join("s1_2.fastq.gz").mksymlinkto(sample1_2)
    p.join("s2_1.fastq.gz").mksymlinkto(sample2_1)
    p.join("s2_2.fastq.gz").mksymlinkto(sample2_2)

    if pytest.config.getoption("-v"):
        logger.info("Created sample flat organization\n\t" + "\n\t".join(str(x) for x in sorted(p.visit())))
    return p
