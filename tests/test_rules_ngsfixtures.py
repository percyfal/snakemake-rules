#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import re
from os.path import basename, realpath, join
import logging
import pytest
import subprocess as sp
from helpers import utils

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

has_ngsfixtures = False

try:
    import pytest_ngsfixtures
    has_ngsfixtures = True
    from pytest_ngsfixtures import factories
    from pytest_ngsfixtures import filetypes
except ImportError as e:
    print("\n\n   pytest-ngsfixtures not installed; install with 'conda install -c percyfal pytest-ngsfixtures'\n\n")

THREADS = pytest.config.getoption("--ngs-threads", "1")
applications = [pytest.config.getoption("--application")] if pytest.config.getoption("--application") else pytest.rules.__all__
#applications = ["samtools"]
rule =  pytest.config.getoption("--rule") if pytest.config.getoption("--rule") else None

def create_testrules(applications, blacklist):
    """Generate list of test rules.

    Params:
      applications (list): list of applications
      blacklist (list): blacklisted rules

    Returns:
      rules (list): list of application directory/rule pairs, where
                    rule is the absolute path to the rule
    """
    testrules = []
    for x in applications:
        for y in getattr(pytest.rules, x):
            if not rule is None:
                if not rule in y:
                    continue
            else:
                if re.sub(".rule", "", basename(y)) in blacklist:
                    continue
            testrules.append((x,realpath(y)))
    return testrules


blacklist_list = []
testrules = create_testrules(applications, blacklist_list)

##############################
# pytest-ngsfixtures
##############################
@pytest.fixture(scope="function", autouse=False)
def data(request, tmpdir_factory):
    """Generate fixture for rule"""
    app, rule = request.param
    if request.function.__name__ == "test_run":
        input, output = utils.parse_rule(rule)
        print(input, output)
        inputfiles = filetypes.filetype_mapping(input)
        fixture_dir = os.path.join(app, os.path.splitext(basename(rule))[0])
        fixture = factories.safe_mktemp(tmpdir_factory, fixture_dir)
        for src in inputfiles:
            p = factories.safe_symlink(fixture, src, basename(src))
    else:
        target = None
        fixture_dir = os.path.join(app, os.path.splitext(basename(rule))[0])
        fixture = factories.safe_mktemp(tmpdir_factory, fixture_dir)
    if request.config.option.ngs_show_fixture:
        logger.info("fixture directory: {}".format(str(fixture)))
        
    return app, rule, target, fixture


@pytest.mark.parametrize("data", sorted(testrules), ids=["{}/{}".format(x[0], basename(x[1])) for x in sorted(testrules)], indirect=["data"])
def test_list(data):
    app, rule, target, fixture = data
    output, err = utils.snakemake_list(fixture, "", snakefile=rule)


blacklist_slow = []
testrules = create_testrules(applications, blacklist_slow)

@pytest.mark.skipif(not has_ngsfixtures, reason="pytest-ngsfixtures not installed; will not run the workflow tests. Install with 'conda install -c percyfal pytest-ngsfixtures'")
@pytest.mark.slow
@pytest.mark.parametrize("data", sorted(testrules), ids=["{}/{}".format(x[0], basename(x[1])) for x in sorted(testrules)], indirect=["data"])
def test_run(data):
    app, rule, target, fixture = data
    if target is None:
        pytest.skip("Unable to parse target for rule {}".format(basename(rule)))
    output, err = utils.snakemake_run(fixture, "", snakefile=rule, target=target)
