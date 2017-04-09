# Copyright (C) 2016 by Per Unneberg
# Helper functions for parsing and for making output executable
import os
from os.path import abspath, dirname, join
import re
import sys
import pytest
import subprocess as sp
import contextlib
import yaml
import logging
from snakemake.parser import parse
from snakemake.io import glob_wildcards, update_wildcard_constraints

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

TESTDIR = join(abspath(dirname(__file__)), os.pardir)
with open(os.path.join(TESTDIR, "rules2target.yaml")) as fh:
    rules2targets = yaml.load(fh)


def run(cmd, stdout=sp.PIPE, stderr=sp.STDOUT):
    """Run workflow subprocess.

    NB: setting env argument to Popen doesn't seem to work; presumably
    the source command wipes env clean

    Args:
      cmd (str): command to run
    """
    close_fds = sys.platform != 'win32'

    proc_prefix = "set -euo pipefail;"
    proc = sp.Popen("{} {} {} {}".format(
        "",
        proc_prefix,
        "",
        cmd),
                    bufsize=-1,
                    shell=True,
                    stdout=stdout,
                    stderr=stderr,
                    close_fds=close_fds,
                    executable="/bin/bash")
    output, err = proc.communicate()
    return output, err

def snakemake_list(fixture, results, **kwargs):
    """Run snakemake list"""
    snakefile = kwargs.get("snakefile", str(fixture.join("Snakefile")))
    args = ['snakemake', '-l', '-d', str(fixture), '-s',
            snakefile]
    save_command(join(str(fixture), "command.sh"), args)
    cmd = " ".join(args)
    output, err = run(cmd, kwargs.get("stdout", sp.PIPE), kwargs.get("stderr", sp.PIPE))
    if not err is None:
        assert(err.decode("utf-8").find("RuleException") == -1), logger.error(err.decode("utf-8"))
        assert(err.decode("utf-8").find("Error:") == -1), logger.error(err.decode("utf-8"))
    if pytest.config.getoption("--show-workflow-output"):
        if not output is None:
            print(output.decode("utf-8"))
        if not err is None:
            print(err.decode("utf-8"))
    return output, err


def snakemake_run(fixture, results, **kwargs):
    """Run snakemake workflow"""
    snakefile = kwargs.get("snakefile", str(fixture.join("Snakefile")))
    targets = kwargs.get("targets", ["all"])
    options = ['-j', kwargs.get("threads", "1"), '-d', str(fixture),
               '-s', snakefile, '--configfile', str(fixture.join('config.yaml'))]

    args = ['snakemake'] + options + targets
    save_command(join(str(fixture), "command.sh"), args)
    cmd = " ".join(args)
    output, err = run(cmd, kwargs.get("stdout", sp.PIPE), kwargs.get("stderr", sp.PIPE))

    if not err is None:
        assert(err.decode("utf-8").find("RuleException") == -1), logger.error(err.decode("utf-8"))
        assert(err.decode("utf-8").find("Error:") == -1), logger.error(err.decode("utf-8"))
    if not output is None:
        print(output.decode("utf-8"))
    if pytest.config.getoption("--show-workflow-output"):
        if not output is None:
            print(output.decode("utf-8"))
        if not err is None:
            print(err.decode("utf-8"))

    # Rerun to get assert statement
    cmd = " ".join(['snakemake'] + options + ['-n'] + targets)
    output, err = run(cmd)
    assert (output.decode("utf-8").find(kwargs.get("results", "Nothing to be done")) >  -1)
    return output, err

# context manager for cd
@contextlib.contextmanager
def cd(path):
    CWD = os.getcwd()
    print("Changing directory from {} to {}".format(CWD, path), file=sys.stderr)

    os.chdir(path)
    try:
        yield
    except:
        print ('Exception caught: ',sys.exc_info()[0], file=sys.stderr)
    finally:
        print("Changing directory back to {}".format(CWD), file=sys.stderr)
        os.chdir(CWD)

def _make_executable(path):
    mode = os.stat(path).st_mode
    mode |= (mode & 0o444) >> 2    # copy R bits to X
    os.chmod(path, mode)

def save_command(fn, args):
    with open(fn, "w") as fh:
        fh.write("#!/bin/bash\n")
        fh.write("PATH={}\n".format(os.environ["PATH"]))
        fh.write("args=$*\n")
        fh.write(" ".join(args) + " ${args}\n")
    _make_executable(fn)


def get_wildcards(inputmap, wildcard_constraints):
    """Given a list of snakemake IO filenames, extract the wildcards.

    Params:
      inputmap (list): list of input wildcard/filename tuples
    """
    d = {}
    for wc, filename in inputmap:
        wc = update_wildcard_constraints(wc, wildcard_constraints, {})
        if filename is None:
            continue
        wildcards = glob_wildcards(wc, [os.path.basename(filename)])
        for k, v in wildcards._asdict().items():
            d[k] = v[0]
    return d


def determine_inputs():
    """Determine the input files.

    Strategy:

    1. try to determine the filetypes, and if successful, see if they
    are mapped to a file
    2. lookup predetermined files in a lookup dictionary
    3. fail, so skip test
    """
    pass
