# -*- coding: utf-8 -*-
import os
import re
from os.path import basename, realpath
import pytest
from helpers import utils, fixture
from snakemake.io import update_wildcard_constraints


def test_wildcards1():
    d = utils.get_wildcards([("{prefix}.bam", "medium.bam")], {})
    assert d['prefix'] == "medium"


def test_wildcards2():
    d = utils.get_wildcards([("{prefix}{ext,.bam}", "medium.bam")], {})
    print(d)

    
def test_determine_fixture():
    # Non-existent filetype
    ft = fixture.determine_fixture("{prefix}.bar")
    assert ft is None
    ft = fixture.determine_fixture("{prefix}.bam")
    assert basename(ft) == "PUR.HG00731.tiny.sort.bam"
    ft = fixture.determine_fixture("config['foo']['bar']")
    assert basename(ft) == "config.yaml"
