# -*- coding: utf-8 -*-
import pytest
from os.path import basename
from helpers import utils, fixture

pytestmark = pytest.mark.skipif(pytest.config.getoption("--application") is not False, reason="application passed; skipping base module tests")


def test_wildcards1():
    d = utils.get_wildcards([('"{prefix}.bam"', "medium.bam")], {})
    assert d['prefix'] == "medium"


def test_wildcards2():
    d = utils.get_wildcards([('"{prefix}{ext,.bam}"', "medium.bam")], {})
    assert d['ext'] == ".bam"


def test_wildcards3():
    d = utils.get_wildcards([('"{prefix}.bar"', "/foo/bar/medium.bar")], {})
    assert d['prefix'] == 'medium'


def test_wildcards4():
    d = utils.get_wildcards([('config[\'foo\'] + ".bar"', "config.yaml")], {})
    assert d == {}


def test_determine_fixture():
    # Non-existent filetype
    ft = fixture.determine_fixture('"{prefix}.bar"')
    assert ft is None
    ft = fixture.determine_fixture('"{prefix}.bam"')
    assert basename(ft) == "PUR.HG00731.tiny.sort.bam"
    ft = fixture.determine_fixture('config[\'foo\'] + ".dict"')
    assert basename(ft) == "scaffolds.dict"
