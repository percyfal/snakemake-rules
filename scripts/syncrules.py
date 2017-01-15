#!/usr/bin/env python3
import os
import sys
import shutil
import argparse
import logging
import difflib
import filecmp
from snakemake_rules import SNAKEMAKE_RULES_PATH
from snakemake.workflow import Workflow
from snakemake.exceptions import print_exception

global snakemake_rule_dict

FORMAT = '%(levelname)s: %(asctime)-15s: %(message)s'
logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(__file__)

def sync_file(source, dest, dryrun=False, diff=False):
    """Sync file source to dest"""
    if diff:
        if not os.path.exists(dest):
            logger.info("Destination '{}' does not exist: skipping diff".format(dest))
            return
        with open(source) as a:
            with open(dest) as b:
                s1 = a.readlines()
                s2 = b.readlines()
                sys.stdout.writelines(difflib.unified_diff(s1, s2, fromfile=source, tofile=dest))
        return
    if not os.path.exists(dest):
        if dryrun:
            logger.info("DRY_RUN: Copying rule '{}' to '{}'".format(source, dest))
        else:
            if not os.path.exists(os.path.dirname(dest)):
                os.makedirs(os.path.dirname(dest))
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
filters = ('.rules', '.rule', '.settings')
snakemake_rule_dict = {}
for path, dirs, files in os.walk(path):
    for f in files:
        if f.endswith(filters):
            if f.endswith(".settings"):
                snakemake_rule_dict[os.path.basename(f)] = os.path.join(path, f)
            else:
                snakemake_rule_dict[os.path.splitext(os.path.basename(f))[0]] = os.path.join(path, f)


parser = argparse.ArgumentParser("Copy/sync rules to a given directory")
parser.add_argument('Snakefile', help="Snakefile to import")
parser.add_argument('-n', '--dry-run', action="store_true", help="Dry run")
parser.add_argument('-d', '--outdir', action="store", default=os.curdir,
                    help="Snakefile to import")
parser.add_argument('-r', '--rule', action="store", default=None, help="rule to sync")
parser.add_argument('-D', '--diff', action="store_true", default=False, help="do diff only")
args = parser.parse_args()

snakefile = os.path.abspath(args.Snakefile)
workflow = create_workflow(snakefile)

# Start work
DEST=args.outdir
workflow_rules = {r.name:r for r in workflow.rules}

for f in [x for x in workflow.included if x.endswith("settings")]:
    basename = os.path.basename(f)
    if args.rule:
        if basename != args.rule:
            continue
    source = snakemake_rule_dict[basename]
    sync_file(source, f, args.dry_run, args.diff)

for r in workflow_rules:
    if args.rule:
        if r != args.rule:
            continue
    if r in list(snakemake_rule_dict.keys()):
        dest = os.path.join(DEST, os.path.relpath(snakemake_rule_dict[r], SNAKEMAKE_RULES_PATH))
        sync_file(snakemake_rule_dict[r], dest, args.dry_run, args.diff)
        rule = create_workflow(snakemake_rule_dict[r]).get_rule(r)
        if rule.conda_env:
            dest = os.path.join(DEST, os.path.relpath(rule.conda_env, SNAKEMAKE_RULES_PATH))
            sync_file(rule.conda_env, dest, args.dry_run, args.diff)
        # Sync also settings file if it exists

    else:
        logger.warn("No such rule '{}' in snakemake_rule_dict".format(r))
