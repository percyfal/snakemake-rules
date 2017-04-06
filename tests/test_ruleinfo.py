# -*- coding: utf-8 -*-
import os
import re
from os.path import basename, realpath
import pytest
from snakemake_rules.core.ruleinfo import _determine_filetypes, parse_rule

    

@pytest.fixture(scope="function")
def rule1(tmpdir_factory):
    rule = """
rule foo:
    input: "{prefix}.foo"
    output: "{prefix}.bar"
    shell: "foo bar {input} {output}"
"""
    p = tmpdir_factory.mktemp("rule1")
    p = p.join("rule1.rule")
    p.write(rule)
    return p
        
def test_parse_rule1(rule1):
    input, output = parse_rule(str(rule1))
    assert input == ['foo']
    assert output == ['"{prefix}.bar"']


@pytest.fixture(scope="function")
def rule2(tmpdir_factory):
    rule = """
rule foo:
    input: foo="{prefix}.foo", bar="{prefix}.bar"
    output: foobar="{prefix}.foobar", barfoo="{prefix}.barfoo"
    shell: "foo bar {input} {output}"
"""
    p = tmpdir_factory.mktemp("rule2")
    p = p.join("rule2.rule")
    p.write(rule)
    return p


def test_parse_rule2(rule2):
    input, output = parse_rule(str(rule2))
    assert sorted(input) == ['bar', 'foo']
    assert sorted(output) == sorted(['{prefix}.foobar', '{prefix}.barfoo'])
    

@pytest.fixture(scope="function")
def rule3(tmpdir_factory):
    rule = """
rule foo:
    input: foo="{prefix}.foo", bar="{prefix}.bar"
    shell: "foo bar {input} {output}"
"""
    p = tmpdir_factory.mktemp("rule3")
    p = p.join("rule3.rule")
    p.write(rule)
    return p


def test_parse_rule3(rule3):
    input, output = parse_rule(str(rule3))
    assert sorted(input) == ['bar', 'foo']
    assert output == []

    
def test_determine_filetypes():
    m = _determine_filetypes('"{prefix}.foo"')
    assert m == ['foo']
    m = _determine_filetypes('foo = "{prefix}.foo", bar = "{prefix}.bar"')
    assert sorted(m) == ['bar', 'foo']
    
