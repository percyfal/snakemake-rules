# Copyright (C) 2016 by Per Unneberg
from os.path import join
import logging
import subprocess as sp
import pytest
from helpers import utils

logger = logging.getLogger(__name__)

stderr = sp.STDOUT
if pytest.config.getoption("--show-workflow-output"):
    stderr = None
THREADS = pytest.config.getoption("--ngs-threads", "1")
skip_workflow = any([pytest.config.getoption("--application"),
                     pytest.config.getoption("--rule")])


@pytest.mark.skipif(
    skip_workflow,
    reason="skipping workflow since --application or --rule passed")
def test_workflow1(snakefile_data):
    wd = str(snakefile_data)
    snakefile = join(wd, 'Snakefile1')
    config = join(wd, 'config1.yaml')
    args = ['snakemake', '-s', snakefile, '-j', str(THREADS), '-d',
            wd, '--configfile', config, "s1.sort.bam.bai"]
    utils.save_command(join(str(wd), "command.sh"), args)
    output = sp.check_output(args, stderr=stderr)


@pytest.mark.skipif(
    skip_workflow,
    reason="skipping workflow since --application or --rule passed")
def test_workflow2(snakefile_data):
    wd = str(snakefile_data)
    snakefile = join(wd, 'Snakefile2')
    config = join(wd, 'config2.yaml')
    args = ['snakemake', '-s', snakefile, '-j', str(THREADS), '-d',
            wd, '--configfile', config]
    utils.save_command(join(str(wd), "command.sh"), args)
    output = sp.check_output(args, stderr=stderr)


# loop different
def test_samples_include(snakefile_data):
    wd = str(snakefile_data)
    snakefile = join(wd, 'Snakefile1')
    config = join(wd, 'config1.yaml')
