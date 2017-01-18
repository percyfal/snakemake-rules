# Copyright (C) 2016 by Per Unneberg
import re
import os
from os.path import abspath, dirname, join, basename
import logging
import shutil
import subprocess as sp
import pytest

logger = logging.getLogger(__name__)

stderr = None if pytest.config.getoption("--show-workflow-output") else sp.STDOUT
applications = [pytest.config.getoption("--application")] if pytest.config.getoption("--application") else pytest.rules.__all__
rule =  pytest.config.getoption("--rule") if pytest.config.getoption("--rule") else None
THREADS = pytest.config.getoption("--threads")


if not set(applications).issubset(pytest.rules.__all__):
    raise Exception("No such application '{}'".format(applications[0]))

blacklist = [
    'picard_do_qc',
    'rsem_calculate_expression_bowtie',
    'snpeff_annotate_variants',
    'snpeff_download_database',
]

rules = []
for x in applications:
    for y in getattr(pytest.rules, x):
        if not rule is None:
            if not rule in y:
                continue
        else:
            if re.sub(".rule", "", basename(y)) in blacklist:
                continue
        rules.append((x,y))


@pytest.mark.parametrize("x", sorted(rules), ids=["{}/{}".format(x[0], basename(x[1])) for x in sorted(rules)])
def test_snakemake_list(x):
    app, rule = x
    name = re.sub(".rule$", "", basename(rule))
    if set([name]).issubset(blacklist):
        pytest.skip("{} part of blacklist".format(name))
    output = sp.check_output(['snakemake', '-s', rule, '-l'], stderr=sp.STDOUT)
    if pytest.config.getoption("--show-workflow-output"):
        print(output.decode("utf-8"))


application_blacklist = ['annovar', 'cloudbiolinux', 'danpos',
                         'dfilter', 'diamond', 'ercc', 'gem', 'homer',
                         'plink', 'snpeff', 'tuxedo', 'utils', 'vcf']
applications = list(set(applications).difference(application_blacklist))

blacklist_slow = [
    'bwa_link_ref',
    'emacs_org_to_reveal',
    'freebayes_parallel',
    'gatk_read_backed_phasing',
    'picard_do_qc',
    'picard_merge_sam',
    'rsem_calculate_expression',
    'rsem_calculate_expression_bowtie',
    'rseqc_clipping_profile',
    'rseqc_qc_8',
    'rseqc_qc',
    'ucsc_download_2bit',
    'ucsc_pseudo',
    'ucsc_link',
    'ucsc_no_alt_analysis_set_reference',
    'ucsc_write_chromosome',
]

slow_rules = []
for x in applications:
    for y in getattr(pytest.rules, x):
        if not rule is None:
            if not rule in y:
                continue
        else:
            if re.sub(".rule", "", basename(y)) in blacklist_slow:
                continue
        slow_rules.append((x,y))


# Helper function to make output executable
def make_executable(path):
    mode = os.stat(path).st_mode
    mode |= (mode & 0o444) >> 2    # copy R bits to X
    os.chmod(path, mode)

@pytest.mark.skipif(not applications, reason="application '{}' in blacklist".format(pytest.config.getoption("--application")))
@pytest.mark.slow
@pytest.mark.parametrize("x", sorted(slow_rules), ids=["{}/{}".format(x[0], basename(x[1])) for x in sorted(slow_rules)])
def test_snakemake_run(x, data):
    app, rule = x
    target = pytest.make_output(rule)
    if target is None:
        pytest.skip("Unable to parse target for rule {}".format(basename(rule)))
    args = ['snakemake', '-f', '-s', rule, '-j', THREADS, '-d', str(data), '--configfile', join(str(data), 'config.yaml')]
    if not target == "config":
        args = args + [target]
    cmdfile = join(str(data), "command.sh")
    with open(cmdfile, "w") as fh:
        fh.write("#!/bin/bash\n")
        fh.write("PATH={}\n".format(os.environ["PATH"]))
        fh.write("args=$*\n")
        fh.write(" ".join(args) + " ${args}\n")
    make_executable(cmdfile)

    output = sp.check_output(args, stderr=stderr)
