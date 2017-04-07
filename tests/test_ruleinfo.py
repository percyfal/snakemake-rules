# -*- coding: utf-8 -*-
import os
import re
from os.path import basename, realpath
import pytest
from snakemake_rules.core.ruleinfo import determine_filetype, parse_rule

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
    d = parse_rule(str(rule1))
    assert d['input'] == ['{prefix}.foo']
    assert d['output'] == ['{prefix}.bar']


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
    d = parse_rule(str(rule2))
    assert sorted(d['input']) == ['{prefix}.bar', '{prefix}.foo']
    assert sorted(d['output']) == sorted(['{prefix}.foobar', '{prefix}.barfoo'])
    

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
    d = parse_rule(str(rule3))
    assert sorted(d['input']) == ['{prefix}.bar', '{prefix}.foo']
    assert d['output'] == []



@pytest.fixture(scope="function")
def rule4(tmpdir_factory):
    rule = """
rule foo:
    input: foo="{prefix}{ext}"
    output: bar="{prefix}(.foo|.bar)"
    shell: "foo bar {input} {output}"
"""
    p = tmpdir_factory.mktemp("rule4")
    p = p.join("rule4.rule")
    p.write(rule)
    return p

def test_parse_rule4(rule4):
    d = parse_rule(str(rule4))
    assert d['input'] == ['{prefix}{ext}']
    assert d['output'] == ['{prefix}(.foo|.bar)']



@pytest.fixture(scope="function")
def rule5(tmpdir_factory):
    rule = """
rule foo:
    input: foo="{prefix}{ext}",
           foobar="{prefix}.foobar"
    output: bar="{prefix}(.foo|.bar)", foo="{prefix}.foo",
            barfoo = "{prefix}.barfoo"
    shell: "foo bar {input} {output}"
"""
    p = tmpdir_factory.mktemp("rule5")
    p = p.join("rule5.rule")
    p.write(rule)
    return p

def test_parse_rule5(rule5):
    d = parse_rule(str(rule5))
    assert sorted(d['input']) == sorted(['{prefix}{ext}', '{prefix}.foobar'])
    assert sorted(d['output']) == sorted(['{prefix}(.foo|.bar)', '{prefix}.foo', '{prefix}.barfoo'])
    
@pytest.fixture(scope="function")
def rule6(tmpdir_factory):
    rule = """
rule foo:
    wildcard_constraints:
      ext = "(.foo|.bar)",
      prefix = "(foobar|barfoo)"
    input: foo="{prefix}{ext}",
           foobar="{prefix}.foobar"
    output: bar="{prefix}{ext}", foo="{prefix}.foo",
            barfoo = "{prefix}.barfoo"
    shell: "foo bar {input} {output}"
"""
    p = tmpdir_factory.mktemp("rule6")
    p = p.join("rule6.rule")
    p.write(rule)
    return p


def test_parse_rule6(rule6):
    d = parse_rule(str(rule6))
    assert sorted(d['input']) == sorted(['{prefix}{ext}', '{prefix}.foobar'])
    assert sorted(d['output']) == sorted(['{prefix}{ext}', '{prefix}.foo', '{prefix}.barfoo'])
    assert d['wildcard_constraints'] == {'ext': '(.foo|.bar)', 'prefix': '(foobar|barfoo)'}
    

    
def test_determine_filetype():
    m = determine_filetype('{prefix}.foo')
    assert m == 'foo'


