# Copyright (C) 2016 by Per Unneberg
# Helper functions for parsning and for making output executable
import os
from os.path import abspath, dirname, join
import re
import sys
import pytest
import subprocess as sp
import contextlib
import yaml
from snakemake.parser import parse

TESTDIR = join(abspath(dirname(__file__)), os.pardir)
with open(os.path.join(TESTDIR, "rules2target.yaml")) as fh:
    rules2targets = yaml.load(fh)

# Map file extension to file extension fixture; keys are concatenated
# with | and compiled to regular expression
re_filetypes = {
    "bam" : "bam",
    "bai" : "bai",
    "bed" : "bed",
    "fa" : "fasta",
    "fasta" : "fasta",
    "fai" : "fastaindex",
    "fofn" : "fofn",
    "sam" : "sam",
    "samtools_stats.txt" : "samtools_stats",
    "vcf" : "vcf",
    "vcf.gz" : "vcf",
    "tbi" : "tabix",
}
regex = re.compile("\.({})\"".format("|".join("{}".format(x) for x in sorted(re_filetypes.keys()))))
regex_input = re.compile("@workflow\.input\(\s+(?P<input>.*)")
regex_output = re.compile("@workflow.output\(\s+(?P<output>.*)")

def parse_rule(rule, prefix=None):
    """Generate input/output information for rule.

    Params:
      rule (str): file containing a snakemake rule """
    rn = os.path.basename(rule).replace(".rule", "")
    app = os.path.basename(os.path.dirname(rule))
    target = rules2targets.get(app, {}).get(rn, None)
    if target:
        return target
    code, linemake, rulecount = parse(rule)
    m = regex_input.search(code)
    if m:
        input = m.group("input")
    m = regex_output.search(code)
    if m is None:
        # return input case
        m = regex_input.search(code)
        output = m.group("input")
    else:
        output = m.group("output")
    m = regex.findall(input)
    if m:
        print(m)
        #print(minput.groups())
    else:
        print("no match")
    m = re.search("\"[ ]*(?P<prefix>\{[a-zA-Z_0-9]+\})+(?P<ext>[_\/\.a-zA-Z0-9 ]+)\"", output)
    # Regular extension; use first one
    if m:
        return "{prefix}{ext}".format(prefix=prefix, ext=m.group("ext"))
    # expand case; skip for now
    m = re.search("expand", output)
    if m:
        return None
    # Config case
    m = re.search("[a-zA-Z =]*(?P<config>config[^\)]+)", output)
    if m:
        return "config"
    return None


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

def snakemake_list(snakemakedata, results, snakefile=None, **kwargs):
    """Run snakemake list"""
    
    d = snakemakedata
    if snakefile is None:
        snakefile = str(d.join("Snakefile"))
    cmd = " ".join(['snakemake', '-l', '-d', str(d), '-s',
                    snakefile])
    save_command(join(str(fixture), "command.sh"), cmd)
    output, err = run(cmd, sp.PIPE, sp.PIPE)
    if not err is None:
        assert(err.decode("utf-8").find("RuleException") == -1), print(err.decode("utf-8"))
    return output, err


def snakemake_run(d, results, snakefile=None, **kwargs):
    """Run snakemake workflow"""
    if snakefile is None:
        snakefile = str(d.join("Snakefile"))
    options = ['-j', kwargs.get("threads", "1"), '-d', str(d),
               '-s', snakefile]
    save_command(join(str(fixture), "command.sh"), cmd)
    output, err = run(cmd, sp.PIPE, sp.PIPE)
    # if kwargs.get("use_conda", False):
    #     options = options + ["--use-conda"]
    cmd = " ".join(['snakemake'] + options + [kwargs.get('target', 'all')])
    cmdfile = os.path.join(str(d), "command.sh")
    with open(cmdfile, "w") as fh:
        fh.write("#!/bin/bash\n")
        fh.write("PATH={}\n".format(os.environ["PATH"]))
        fh.write("args=$*\n")
        fh.write(cmd + " ${args}\n")
    make_executable(cmdfile)

    output, err = run(cmd, kwargs.get("stdout", None), kwargs.get("stderr", None))
    # Rerun to get assert statement
    cmd = " ".join(['snakemake'] + options + ['-n'] + [kwargs.get('target', 'all')])
    output, err = run(cmd)
    assert (output.decode("utf-8").find(kwargs.get("results", "Nothing to be done")) >  -1)
    

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

