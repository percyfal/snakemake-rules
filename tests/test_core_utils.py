# -*- coding: utf-8 -*-
import os
import re
from os.path import basename, realpath
import pytest
from snakemake_rules.core import utils
from snakemake.logging import logger


@pytest.fixture(autouse=False)
def sampleinfo(tmpdir_factory):
    p = tmpdir_factory.mktemp("config")
    p = p.join("sampleinfo.csv")
    s = """SM,PU,fastq
s1,pu1,fastq1
s2,pu1,fastq1
s2,pu2,fastq2
s3,pu1,fastq3
"""
    p.write(s)
    return p


@pytest.fixture(scope="function", autouse=False,
                params=[None, "sampleinfo", "wrong sampleinfo", "no include"])
def config(request, sampleinfo):
    config = {'samples': ['s1', 's2'], 'ignore_samples': ['s3'], 'settings': {}}
    if request.param in ["sampleinfo", "no include"]:
        config['settings']['sampleinfo'] = str(sampleinfo)
    if request.param == "wrong sampleinfo":
        config['settings']['sampleinfo'] = "foo.csv"
    if request.param == "no include":
        config['samples'] = []
    return config, request.param


@pytest.mark.skipif(pytest.config.getoption("--application") is not False, reason="application passed; skipping core utils tests")
def test_config(config):
    conf, param = config
    utils.get_samples(conf, logger)
    if param is None:
        assert conf.get('_sampleinfo') is None
        assert conf['settings'].get('sampleinfo') is None
    elif param in ["sampleinfo", "no include"]:
        assert conf.get('_sampleinfo') is not None
        assert len(conf['_sampleinfo']) == 3
        assert len(conf['samples']) == 2
    else:
        assert conf['settings']['sampleinfo'] == 'foo.csv'
