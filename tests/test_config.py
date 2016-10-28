# Copyright (C) 2016 by Per Unneberg
import os
from os.path import join
import shutil
import subprocess as sp
import pytest

def test_conditional_include():
    output = sp.check_output(['snakemake', '-s', join(pytest.TESTDIR, os.pardir, "bwa", "bwa_index.rule"), '-l'], stderr=sp.STDOUT)
    print(output.decode("utf-8"))
