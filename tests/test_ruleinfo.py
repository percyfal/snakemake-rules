# -*- coding: utf-8 -*-
import os
import re
from os.path import basename, realpath
import pytest
from snakemake_rules.core.ruleinfo import parse_rule, regex_target_list, get_targets


pytestmark = pytest.mark.skipif(pytest.config.getoption("--application") is not False, reason="application passed; skipping ruleinfo tests")


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
    assert d['input']['_list'] == ['"{prefix}.foo"']
    assert d['output']['_list'] == ['"{prefix}.bar"']


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
    assert d['input'] == {'bar': '"{prefix}.bar"', 'foo': '"{prefix}.foo"', '_list': []}
    assert d['output'] == {'foobar': '"{prefix}.foobar"', 'barfoo': '"{prefix}.barfoo"', '_list': []}


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
    assert d['input'] == {'bar': '"{prefix}.bar"', 'foo': '"{prefix}.foo"', '_list': []}
    assert d['output'] == {'_list': []}



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
    assert d['input'] == {'foo': '"{prefix}{ext}"', '_list': []}
    assert d['output'] == {'bar': '"{prefix}(.foo|.bar)"', '_list': []}



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
    assert d['input'] == {'foo': '"{prefix}{ext}"', 'foobar': '"{prefix}.foobar"', '_list': []}
    assert d['output'] == {'bar': '"{prefix}(.foo|.bar)"', 'foo': '"{prefix}.foo"',
                           'barfoo': '"{prefix}.barfoo"', '_list': []}


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
    assert d['input'] ==  {'foo': '"{prefix}{ext}"', 'foobar': '"{prefix}.foobar"', '_list': []}
    assert d['output'] == {'bar': '"{prefix}{ext}"', 'foo': '"{prefix}.foo"',
                           'barfoo': '"{prefix}.barfoo"', '_list': []}
    assert d['wildcard_constraints'] == {'ext': '(.foo|.bar)', 'prefix': '(foobar|barfoo)'}


@pytest.fixture(scope="function")
def rule7(tmpdir_factory):
    rule = """
rule foo:
    wildcard_constraints:
      ext = "(.foo|.bar)",
      prefix = "(foobar|barfoo)"
    input: foo="{prefix}{ext}",
           foobar="{prefix}.foobar",
           ref = config['foo']['ref'],
           fai = config['foo']['ref'] + ".fai",
    output: bar="{prefix}{ext}", foo="{prefix}.foo",
            barfoo = "{prefix}.barfoo"
    shell: "foo bar {input} {output}"
"""
    p = tmpdir_factory.mktemp("rule7")
    p = p.join("rule7.rule")
    p.write(rule)
    return p


def test_parse_rule7(rule7):
    d = parse_rule(str(rule7))
    assert d['input'] ==  {'foo': '"{prefix}{ext}"', 'foobar': '"{prefix}.foobar"', '_list': [],
                           'ref': 'config[\'foo\'][\'ref\']', 'fai': 'config[\'foo\'][\'ref\']+".fai"'}


@pytest.fixture(scope="function")
def rule8(tmpdir_factory):
    rule = """
rule foo:
    input: "{prefix}.bar", config['foo']['bar']
    output: foo="{prefix}.foo",
    shell: "foo bar {input} {output}"
"""
    p = tmpdir_factory.mktemp("rule8")
    p = p.join("rule8.rule")
    p.write(rule)
    return p


def test_parse_rule8(rule8):
    d = parse_rule(str(rule8))
    assert d['input'] ==  {'_list': ['"{prefix}.bar"', 'config[\'foo\'][\'bar\']']}



def test_get_targets():
    d = get_targets("\"{prefix}.bar\",foo=\"{prefix}.foo\"")
    assert d['_list'] == ['"{prefix}.bar"']
    assert d['foo'] == '"{prefix}.foo"'
    d = get_targets("foo=config['foo']['bar']")
    assert d['foo'] == "config['foo']['bar']"
    d = get_targets("foo=os.path.splitext(config['foo']['bar'])[0]+\".bar\"")
    assert d['foo'] == 'os.path.splitext(config[\'foo\'][\'bar\'])[0]+".bar"'
