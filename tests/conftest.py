# Copyright (C) 2016 by Per Unneberg
import os
from os.path import abspath, dirname, join, isdir
import sys
import re
import logging
import pytest
import shutil
import subprocess as sp
import ast
import yaml
from snakemake.parser import parse

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

TESTDIR = abspath(dirname(__file__))
RULEDIR = join(abspath(dirname(__file__)), os.pardir, "snakemake_rules")
SNAKEFILE = join(TESTDIR, "Snakefile")
SNAKEFILE_REGIONS = join(TESTDIR, "Snakefile_regions")
CONFIG = join(TESTDIR, "config.yaml")
CONFIG_REGIONS = join(TESTDIR, "config_regions.yaml")

# Add test source path to pythonpath
sys.path.insert(0, join(TESTDIR, os.pardir))

# Autogenerate all application directories
blacklist = ['__pycache__', 'bio', 'comp']
applications = {x:join(RULEDIR, x) for x in os.listdir(RULEDIR) if isdir(join(RULEDIR, x)) and not x in blacklist}
rules = {k: [join(RULEDIR, k, x) for x in os.listdir(v) if x.endswith(".rule")] for k,v in applications.items() }
with open(os.path.join(TESTDIR, "rules2target.yaml")) as fh:
    rules2targets = yaml.load(fh)


def make_output(rule, prefix="s1"):
    rn = os.path.basename(rule).replace(".rule", "")
    app = os.path.basename(os.path.dirname(rule))
    target = rules2targets.get(app, {}).get(rn, None)
    if target:
        return target
    code, linemake, rulecount = parse(rule)
    m = re.search("@workflow.output\(\s+(?P<output>.*)", code)
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
                     help="application to test",
                     dest="application")
    parser.addoption("-T", "--threads", action="store",
                     default="1",
                     help="number of threads to use",
                     dest="threads")
    parser.addoption("-R", "--rule", action="store",
                     default=None,
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
        'snakefile': SNAKEFILE,
        'snakefile_regions' : SNAKEFILE_REGIONS,
        'config' : CONFIG,
        'config_regions' : CONFIG_REGIONS,
        'rules' : rules,
        'make_output' : make_output,
    }

##################################################
# Setup test fixtures
##################################################
# Input files

# metadata
sampleinfo = abspath(join(dirname(__file__), "data", "sampleinfo.csv"))
configfile = abspath(join(dirname(__file__), "data", "config.yaml"))

# references
chr11 = abspath(join(dirname(__file__), "data", "chr11.fa"))
chr11fai = abspath(join(dirname(__file__), "data", "chr11.fa.fai"))
chr11dict = abspath(join(dirname(__file__), "data", "chr11.dict"))
chromsizes = abspath(join(dirname(__file__), "data", "chrom.sizes"))

# annotation files
dbsnp = abspath(join(dirname(__file__), "data", "dbsnp132_chr11.vcf"))
ref_transcripts = abspath(join(dirname(__file__), "data", "ref-transcripts.gtf"))
ref_transcripts_bed12 = abspath(join(dirname(__file__), "data", "ref-transcripts.bed12"))
ref_transcripts_genepred = abspath(join(dirname(__file__), "data", "ref-transcripts.genePred"))
targets = abspath(join(dirname(__file__), "data", "targets.bed"))
targets_list = abspath(join(dirname(__file__), "data", "targets.interval_list"))

# fastq files
sample1_1 = abspath(join(dirname(__file__), "data", "s1_1.fastq.gz"))
sample1_2 = abspath(join(dirname(__file__), "data", "s1_2.fastq.gz"))
sample2_1 = abspath(join(dirname(__file__), "data", "s2_1.fastq.gz"))
sample2_2 = abspath(join(dirname(__file__), "data", "s2_2.fastq.gz"))


# alignment files
sam = abspath(join(dirname(__file__), "data", "s1.sam"))
bam = abspath(join(dirname(__file__), "data", "s1.bam"))
sortbam = abspath(join(dirname(__file__), "data", "s1.sort.bam"))
sortbed = abspath(join(dirname(__file__), "data", "s1.sort.bed"))
sortwig = abspath(join(dirname(__file__), "data", "s1.sort.wig"))
sortbedgraph = abspath(join(dirname(__file__), "data", "s1.sort.bedGraph"))
sortbambai = abspath(join(dirname(__file__), "data", "s1.sort.bam.bai"))
rgbam = abspath(join(dirname(__file__), "data", "s1.rg.bam"))
rgsortbam = abspath(join(dirname(__file__), "data", "s1.rg.sort.bam"))
s2rgsortbam = abspath(join(dirname(__file__), "data", "s2.rg.sort.bam"))
bamfofn = abspath(join(dirname(__file__), "data", "bamfiles.fofn"))

# vcf
vcf = abspath(join(dirname(__file__), "data", "s1.vcf"))
vcfgz = abspath(join(dirname(__file__), "data", "s1.vcf.gz"))
vcffofn = abspath(join(dirname(__file__), "data", "s1.vcf.fofn"))
gvcf = abspath(join(dirname(__file__), "data", "s1.g.vcf"))
gvcfgz = abspath(join(dirname(__file__), "data", "s1.g.vcf.gz"))


##############################
# sample
##############################
@pytest.fixture(scope="function", autouse=False)
def data(tmpdir_factory):
    """
    Setup input data
    """
    p = tmpdir_factory.mktemp('data')

    # metadata
    p.join("sampleinfo.csv").mksymlinkto(sampleinfo)
    p.join("config.yaml").mksymlinkto(configfile)

    # references
    p.join("chr11.fa").mksymlinkto(chr11)
    p.join("chr11.fa.fai").mksymlinkto(chr11fai)
    p.join("chr11.dict").mksymlinkto(chr11dict)
    p.join("chrom.sizes").mksymlinkto(chromsizes)

    # annotation files
    p.join("dbsnp132_chr11.vcf").mksymlinkto(dbsnp)
    p.join("ref-transcripts.gtf").mksymlinkto(ref_transcripts)
    p.join("s1.gtf").mksymlinkto(ref_transcripts)
    p.join("ref-transcripts.bed12").mksymlinkto(ref_transcripts_bed12)
    p.join("ref-transcripts.genePred").mksymlinkto(ref_transcripts_genepred)
    p.join("s1.genePred").mksymlinkto(ref_transcripts_genepred)
    p.join("targets.bed").mksymlinkto(targets)
    p.join("targets.interval_list").mksymlinkto(targets_list)

    # fasta files
    p.join("s1.fasta").mksymlinkto(chr11)
    p.join("s1.fa").mksymlinkto(chr11)
    
    # fastq files
    p.join("s1.fastq.gz").mksymlinkto(sample1_1)
    p.join("s1_1.fastq.gz").mksymlinkto(sample1_1)
    p.join("s1_2.fastq.gz").mksymlinkto(sample1_2)
    p.join("s2_1.fastq.gz").mksymlinkto(sample2_1)
    p.join("s2_2.fastq.gz").mksymlinkto(sample2_2)

    # alignment files
    p.join("s1.sam").mksymlinkto(sam)
    p.join("s1.bam").mksymlinkto(bam)
    p.join("s1.bed").mksymlinkto(sortbed)
    p.join("s1.sort.bam").mksymlinkto(sortbam)
    p.join("s2.sort.bam").mksymlinkto(sortbam)
    p.join("s1.sort.bam.bai").mksymlinkto(sortbambai)
    p.join("s2.sort.bam.bai").mksymlinkto(sortbambai)
    p.join("s1.bdg").mksymlinkto(sortbedgraph)
    p.join("s1.wig").mksymlinkto(sortwig)
    p.join("s1.rg.bam").mksymlinkto(rgbam)
    p.join("s1.rg.sort.bam").mksymlinkto(rgsortbam)
    p.join("s2.rg.sort.bam").mksymlinkto(s2rgsortbam)
    p.join("bamfiles.fofn").mksymlinkto(bamfofn)
    p.join("s1.bam.fofn").mksymlinkto(bamfofn)

    # vcf files
    p.join("s1.vcf").mksymlinkto(vcf)
    p.join("s1.g.vcf").mksymlinkto(gvcf)
    p.join("s1.g.vcf.gz").mksymlinkto(gvcfgz)
    p.join("s1.vcf.fofn").mksymlinkto(vcffofn)
    p.join("s1.g.vcf.fofn").mksymlinkto(vcffofn)
    p.join("s2.g.vcf").mksymlinkto(gvcf)
    p.join("s2.g.vcf.gz").mksymlinkto(gvcfgz)

    return p
