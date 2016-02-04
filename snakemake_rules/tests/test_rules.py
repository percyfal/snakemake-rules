# Copyright (C) 2016 by Per Unneberg
from os.path import abspath, dirname, join
import logging
import shutil
import subprocess as sp
import pytest

logger = logging.getLogger(__name__)

SNAKEFILE=join(abspath(dirname(__file__)), "Snakefile")
SNAKEFILE_REGIONS=join(abspath(dirname(__file__)), "Snakefile_regions")


def test_snakemake():
    """Test snakemake command call"""
    output = sp.check_output(['snakemake', '-s', SNAKEFILE, '-l'], stderr=sp.STDOUT)
    assert "bwa_mem" in output.decode("utf-8")
    

def test_bwa_align():
    """Test bwa alignment"""
    output = sp.check_output(['snakemake', '-s', SNAKEFILE, '-F', 'data/test.sort.bam'], stderr=sp.STDOUT)
    assert "Removing temporary output file data/chr11.fa.sa." in output.decode("utf-8").replace("\t", " ")
    assert "Removing temporary output file data/test.bam." in output.decode("utf-8").replace("\t", " ")
    assert "3 of 3 steps (100%) done" in output.decode("utf-8").replace("\t", " ")



def test_bamtools_filter():
    """Test bamtools filter without using script file, dry run"""
    output = sp.check_output(['snakemake', '-s', SNAKEFILE, '-F',  '-n',  '-p', 'data/test.filter.bam'])
    assert 'bamtools filter -in data/test.bam -out data/test.filter.bam -mapQuality ">=255"   > data/test.filter.log' in output.decode("utf-8").replace("\t", " ")


def test_bamtools_filter_script():
    """Test bamtools filter using script file, dry run. See issue #12."""
    output = sp.check_output(['snakemake', '-s', SNAKEFILE_REGIONS, '-F',  '-n',  '-p', 'data/test.filter.bam'])
    s = "echo '{\"filters\":[\n{\n\"reference\":\"chr11\",\n \"mapQuality\": \">=255\"\n\n}\n    ]\n}\n' > data/test.script"
    assert s in output.decode("utf-8")
    s = 'bamtools filter -in data/test.bam -out data/test.filter.bam -mapQuality ">=255" -script  data/test.script > data/test.filter.log'
    assert s in output.decode("utf-8")
