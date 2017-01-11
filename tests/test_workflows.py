# Copyright (C) 2016 by Per Unneberg
import re
from os.path import abspath, dirname, join, basename
import logging
import shutil
import subprocess as sp
import pytest

logger = logging.getLogger(__name__)

stderr = None if pytest.config.getoption("--show-workflow-output") else sp.STDOUT
THREADS = pytest.config.getoption("--threads")


def test_workflow1(snakefile_data):
    wd = str(snakefile_data)
    snakefile = join(wd, 'Snakefile1')
    config = join(wd, 'config1.yaml')
    args = ['snakemake', '-s', snakefile, '-j', THREADS, '-d',
            wd, '--configfile', config, 's1.sort.bai']
    output = sp.check_output(args, stderr=stderr)

    
def test_workflow2(snakefile_data):
    wd = str(snakefile_data)
    snakefile = join(wd, 'Snakefile2')
    config = join(wd, 'config2.yaml')
    args = ['snakemake', '-s', snakefile, '-j', THREADS, '-d',
            wd, '--configfile', config]#, 'all']
    output = sp.check_output(args, stderr=stderr)
