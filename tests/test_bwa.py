# Copyright (C) 2016 by Per Unneberg
from os.path import join, basename
import shutil
import subprocess as sp
import pytest

stderr = None if pytest.config.getoption("--show-workflow-output") else sp.STDOUT

bwa = pytest.mark.skipif(shutil.which("bwa") is None, reason="bwa not installed")
samtools = pytest.mark.skipif(shutil.which("samtools") is None, reason="samtools not installed")

blacklist = []
rules = [(x) for x in pytest.rules.bwa if not basename(x).rsplit(".rule") in blacklist]

@pytest.fixture(scope="function", autouse=False)
def config(data):
    p = data.join("config.yaml")
    p.write("bwa:\n  ref: chr11.fa\n  index: chr11.fa")


@pytest.mark.parametrize("rule", rules)
def test_bwa_list(rule):
    output = sp.check_output(['snakemake', '-s', rule, '-l'], stderr=sp.STDOUT)
    
@bwa
@samtools
@pytest.mark.slow
@pytest.mark.parametrize("rule", rules)
def test_bwa_run(rule, data, config):
    target = pytest.make_output(rule)
    if target is None:
        pytest.skip("Unable to parse target for rule {}".format(basename(rule)))
    args = ['snakemake', '-s', rule, '-d', str(data), '--configfile', join(str(data), 'config.yaml')]
    if not target == "config":
        args = args + [target]
    output = sp.check_output(args, stderr=stderr)
