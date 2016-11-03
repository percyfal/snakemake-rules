# Copyright (C) 2016 by Per Unneberg
import re
from os.path import abspath, dirname, join, basename
import logging
import shutil
import subprocess as sp
import pytest

logger = logging.getLogger(__name__)

stderr = None if pytest.config.getoption("--show-workflow-output") else sp.STDOUT
applications = [pytest.config.getoption("--application")] if pytest.config.getoption("--application") else pytest.rules.__all__

if not set(applications).issubset(pytest.rules.__all__):
    raise Exception("No such application '{}'".format(applications[0]))


blacklist =['rsem_calculate_expression_bowtie']
rules = [(x, y) for x in applications for y in getattr(pytest.rules, x) if not basename(y).rsplit(".rule") in blacklist]

@pytest.mark.parametrize("x", rules, ids=["{}/{}".format(x[0], basename(x[1])) for x in rules])
def test_snakemake_list(x):
    app, rule = x
    name = re.sub(".rule$", "", basename(rule))
    if set([name]).issubset(blacklist):
        pytest.skip("{} part of blacklist".format(name))
    output = sp.check_output(['snakemake', '-s', rule, '-l'], stderr=sp.STDOUT)


application_blacklist = ['utils']
applications = list(set(applications).difference(application_blacklist))
rules = [(x, y) for x in applications for y in getattr(pytest.rules, x) if not basename(y).rsplit(".rule") in blacklist]

@pytest.mark.skipif(not applications, reason="application '{}' in blacklist".format(pytest.config.getoption("--application")))
@pytest.mark.slow
@pytest.mark.parametrize("x", rules, ids=["{}/{}".format(x[0], basename(x[1])) for x in rules])
def test_snakemake_run(x, data):
    app, rule = x
    target = pytest.make_output(rule)
    if target is None:
        pytest.skip("Unable to parse target for rule {}".format(basename(rule)))
    args = ['snakemake', '-s', rule, '-d', str(data), '--configfile', join(str(data), 'config.yaml')]
    if not target == "config":
        args = args + [target]
    output = sp.check_output(args, stderr=stderr)



