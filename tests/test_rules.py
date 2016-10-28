# Copyright (C) 2016 by Per Unneberg
from os.path import abspath, dirname, join
import logging
import shutil
import subprocess as sp
import pytest

logger = logging.getLogger(__name__)





# SNAKEFILE=join(abspath(dirname(__file__)), "Snakefile")
# SNAKEFILE_REGIONS=join(abspath(dirname(__file__)), "Snakefile_regions")
# DRYRUN="-n"

# def test_snakemake():
#     """Test snakemake command call"""
#     output = sp.check_output(['snakemake', '-s', SNAKEFILE, '-l'], stderr=sp.STDOUT)
#     assert "bwa_mem" in output.decode("utf-8")
    
# @pytest.mark.skipif(shutil.which("bwa") is None, reason="bwa not installed")
# def test_bwa_align():
#     """Test bwa alignment"""
#     output = sp.check_output(['snakemake', '-s', SNAKEFILE, '-F', 'data/s1.sort.bam'], stderr=sp.STDOUT)
#     assert "Removing temporary output file data/chr11.fa.sa." in output.decode("utf-8").replace("\t", " ")
#     assert "Removing temporary output file data/s1.bam." in output.decode("utf-8").replace("\t", " ")
#     assert "3 of 3 steps (100%) done" in output.decode("utf-8").replace("\t", " ")



# def test_bamtools_filter():
#     """Test bamtools filter without using script file, dry run"""
#     output = sp.check_output(['snakemake', '-s', SNAKEFILE, '-F',  '-n',  '-p', 'data/s1.filter.bam'])
#     assert 'bamtools filter -in data/s1.bam -out data/s1.filter.bam -mapQuality ">=255"   > data/s1.filter.log' in output.decode("utf-8").replace("\t", " ")


# def test_bamtools_filter_script():
#     """Test bamtools filter using script file, dry run. See issue #12."""
#     output = sp.check_output(['snakemake', '-s', SNAKEFILE_REGIONS, '-F',  '-n',  '-p', 'data/s1.filter.bam'])
#     s = "echo '{\"filters\":[\n{\n\"reference\":\"chr11\",\n \"mapQuality\": \">=255\"\n\n}\n    ]\n}\n' > data/s1.script"
#     assert s in output.decode("utf-8")
#     s = 'bamtools filter -in data/s1.bam -out data/s1.filter.bam -mapQuality ">=255" -script  data/s1.script > data/s1.filter.log'
#     assert s in output.decode("utf-8")


# def test_picard_merge():
#     """Test picard merge."""
#     output = sp.check_output(['snakemake', '-s', SNAKEFILE, '-F',  '-n',  '-p', 'data/s.merge.bam'])
#     assert 'input: data/s1.bam, data/s2.bam' in output.decode("utf-8")
#     assert 'output: data/s.merge.bam' in output.decode("utf-8")
        
