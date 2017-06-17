#!/usr/bin/env python3
import os
import sys
from os.path import abspath, dirname, join, relpath, exists, splitext
import shutil
import argparse
import logging
import difflib
import filecmp
from snakemake_rules import SNAKEMAKE_RULES_PATH
from snakemake.workflow import Workflow
from snakemake.exceptions import print_exception

FORMAT = '%(levelname)s: %(asctime)-15s: %(message)s'
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__file__)

global snakemake_rule_dict


def sync_file(source, dest, dryrun=False, diff=False):
    """Sync file source to dest"""
    if diff:
        if not exists(dest):
            logger.info("Destination '{}' does not exist: skipping diff".format(dest))
            return
        with open(source) as a:
            with open(dest) as b:
                s1 = a.readlines()
                s2 = b.readlines()
                sys.stdout.writelines(difflib.unified_diff(s1, s2, fromfile=source, tofile=dest))
        return
    if not exists(dest):
        if dryrun:
            logger.info("DRY_RUN: Copying rule '{}' to '{}'".format(source, dest))
        else:
            if not exists(dirname(dest)):
                os.makedirs(dirname(dest))
            logger.info("Copying rule '{}' to '{}'".format(source, dest))
            shutil.copy2(source, dest)
    else:
        equal = filecmp.cmp(source, dest)
        if (not equal):
            if dryrun:
                logger.info("DRY_RUN: Updating rule '{}' to '{}'".format(source, dest))
            else:
                logger.info("Updating rule '{}' to '{}'".format(source, dest))
                shutil.copy2(source, dest)
        else:
            if dryrun:
                logger.info("DRY_RUN: rule '{}' up to date".format(dest))
            else:
                logger.info("rule '{}' up to date".format(dest))


def create_workflow(snakefile):
    workflow = Workflow(snakefile=snakefile, use_conda=True)

    try:
        workflow.include(snakefile,
                         overwrite_first_rule=True,
                         print_compilation=False)
        workflow.check()
    except (Exception, BaseException) as ex:
        print_exception(ex, workflow.linemaps)
        success = False

    return workflow


# Snakemake rules dictionary
path = SNAKEMAKE_RULES_PATH
filters = ('.rules', '.rule', '.settings', '.py')
snakemake_rule_dict = {}
for path, dirs, files in os.walk(path):
    for f in files:
        if f.endswith(filters):
            if f.endswith(".settings"):
                snakemake_rule_dict[f] = join(path, f)
            elif f.endswith(".py"):
                mod = ".".join(relpath(join(path, f), SNAKEMAKE_RULES_PATH).split(os.sep))
                snakemake_rule_dict[mod] = join(path, f)
            else:
                rule = splitext(f)[0]
                snakemake_rule_dict[rule] = join(path, f)

USAGE = """Copy/sync rules to a given directory"""
EPILOG = """

Run syncrules.py on a Snakefile and sync snakemake_rules rules to
directory provided by option -d. If the Snakefile loads configuration
files syncrules.py may not be able to find them if the Snakefile does
not reside in the working directory. Consequently it may be necessary
to change directory to that of the input Snakefile.
"""


parser = argparse.ArgumentParser(USAGE, epilog=EPILOG)
parser.add_argument('Snakefile', help="Snakefile to import")
parser.add_argument('-n', '--dry-run', action="store_true", help="Dry run")
parser.add_argument('-d', '--outdir', action="store", default=os.curdir,
                    help="output directory")
parser.add_argument('-r', '--rule', action="store", default=None, help="rule to sync")
parser.add_argument('-m', '--module', action="store", default=None, help="python module to sync")
parser.add_argument('-s', '--settings', action="store", default=None, help="settings file to sync")
parser.add_argument('-D', '--diff', action="store_true", default=False, help="do diff only")
parser.add_argument('-v', '--verbose', action="store_true", default=False, help="increase verbosity")
args = parser.parse_args()

if args.verbose:
    logger.setLevel(logging.DEBUG)

snakefile = abspath(args.Snakefile)
workflow = create_workflow(snakefile)

# Start work
DEST = args.outdir
workflow_rules = {r.name: r for r in workflow.rules}

# Sync python module
if args.module:
    assert args.module.endswith(".py"), "the supplied module argument '{}' misses the suffix '.py'".format(args.module)
    try:
        source = snakemake_rule_dict[args.module]
    except KeyError as e:
        print("No such module '{}'".format(args.module))
        raise e
    dest = join(DEST, args.module)
    sync_file(source, dest, args.dry_run, args.diff)
    sys.exit(0)

# Sync separate rule
if args.rule:
    try:
        source = snakemake_rule_dict[args.rule]
    except KeyError as e:
        print("No such rule '{}'".format(args.rule))
        raise e
    dest = join(DEST, relpath(source, SNAKEMAKE_RULES_PATH))
    sync_file(source, dest, args.dry_run, args.diff)
    sys.exit(0)

# Sync settings
if args.settings:
    try:
        source = snakemake_rule_dict[args.settings]
    except KeyError as e:
        print("No such rule '{}'".format(args.settings))
        raise e
    dest = join(DEST, relpath(source, SNAKEMAKE_RULES_PATH))
    sync_file(source, dest, args.dry_run, args.diff)
    sys.exit(0)

# Sync all rules in workflow
envfiles = []
for r in workflow_rules:
    if args.rule:
        if r != args.rule:
            continue
    if r in list(snakemake_rule_dict.keys()):
        source = snakemake_rule_dict[r]
        dest = join(DEST, relpath(snakemake_rule_dict[r], SNAKEMAKE_RULES_PATH))
        sync_file(source, dest, args.dry_run, args.diff)
        rule = create_workflow(snakemake_rule_dict[r]).get_rule(r)
        env = join(dirname(snakemake_rule_dict[r]), "env.yaml") if exists(join(dirname(snakemake_rule_dict[r]), "env.yaml")) else rule.conda_env
        if env in envfiles:
            continue
        if env:
            dest = join(DEST, relpath(env, SNAKEMAKE_RULES_PATH))
            sync_file(env, dest, args.dry_run, args.diff)
            envfiles.append(env)

    else:
        logger.warn("No such rule '{}' in snakemake_rule_dict".format(r))
