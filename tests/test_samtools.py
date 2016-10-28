# Copyright (C) 2016 by Per Unneberg
from os.path import join, basename
import shutil
import subprocess as sp
import pytest

stderr = None if pytest.config.getoption("--show-workflow-output") else sp.STDOUT

samtools = pytest.mark.skipif(shutil.which("samtools") is None, reason="samtools not installed")

blacklist = []
rules = [(x) for x in pytest.rules.samtools if not basename(x).rsplit(".rule") in blacklist]

@pytest.mark.parametrize("rule", rules)
def test_samtools_list(rule):
    output = sp.check_output(['snakemake', '-s', rule, '-l'], stderr=sp.STDOUT)

@samtools
@pytest.mark.slow
@pytest.mark.parametrize("rule", rules)
def test_samtools_run(rule, data):
    target = pytest.make_output(rule)
    if target is None:
        pytest.skip("Unable to parse target for rule {}".format(basename(rule)))
    args = ['snakemake', '-s', rule, '-f', '-d', str(data)]
    if not target == "config":
        args = args + [target]
    output = sp.check_output(args, stderr=stderr)
