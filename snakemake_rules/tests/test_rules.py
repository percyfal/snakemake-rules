# Copyright (C) 2016 by Per Unneberg
from os.path import abspath, dirname, join
import logging
import shutil
import subprocess as sp
import pytest

logger = logging.getLogger(__name__)

SNAKEFILE=join(abspath(dirname(__file__)), "Snakefile")


def test_snakemake():
    """Test snakemake command call"""
    output = sp.check_output(['snakemake', '-s', SNAKEFILE, '-l'], stderr=sp.STDOUT)
    assert "bwa_mem" in str(output)
    

def test_bwa_align():
    """Test bwa alignment"""
    output = sp.check_output(['snakemake', '-s', SNAKEFILE, '-F', 'data/test.sort.bam'], stderr=sp.STDOUT)
    assert "3 of 3 steps (100%) done" in str(output)
