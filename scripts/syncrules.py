#!/usr/bin/env python3
import os
import shutil
import argparse
import logging
from snakemake_rules import SNAKEMAKE_RULES_PATH
from snakemake.workflow import Workflow
from snakemake.exceptions import print_exception

FORMAT = '%(levelname)s: %(asctime)-15s: %(message)s'
logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(__file__)

def sync_file(source, dest, dryrun=False):
    """Sync file source to dest"""
    if not os.path.exists(dest):
        if dryrun:
            logger.info("DRY_RUN: Copying rule '{}' to '{}'".format(source, dest))
        else:
            if not os.path.exists(os.path.dirname(dest)):
                os.makedirs(os.path.dirname(dest))
            logger.info("Copying rule '{}' to '{}'".format(source, dest))
            shutil.copy2(source, dest)
    else:
        srctime = os.path.getmtime(source)
        desttime = os.path.getmtime(dest)
        if (desttime > srctime):
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

# Snakemake rules dictionary
path = SNAKEMAKE_RULES_PATH
filters = ('.rules', '.rule', '.settings')
snakemake_rule_dict = {}
for path, dirs, files in os.walk(path):
    for f in files:
        if f.endswith(filters):
            snakemake_rule_dict[os.path.splitext(os.path.basename(f))[0]] = os.path.join(path, f)


parser = argparse.ArgumentParser("Copy/sync rules to a given directory")
parser.add_argument('Snakefile', help="Snakefile to import")
parser.add_argument('-n', '--dry-run', action="store_true", help="Dry run")
parser.add_argument('-d', '--outdir', action="store", default=os.curdir,
                    help="Snakefile to import")
args = parser.parse_args()

snakefile = os.path.abspath(args.Snakefile)

workflow = Workflow(snakefile=snakefile)

try:
    workflow.include(snakefile,
                     overwrite_first_rule=True,
                     print_compilation=False)
    workflow.check()
except (Exception, BaseException) as ex:
    print_exception(ex, workflow.linemaps)
    success = False

# Map the rules included from snakemake_rules
DEST=args.outdir

included_rules = {x:os.path.join(DEST, os.path.relpath(x, SNAKEMAKE_RULES_PATH)) for x in workflow.included if x.startswith(SNAKEMAKE_RULES_PATH)}
# Copy rules to outdir
for k, v in included_rules.items():
    sync_file(k, v, args.dry_run)

workflow_rules = {r.name:r for r in workflow.rules}

leftover_rules = set(workflow_rules.keys()).difference(included_rules.keys())
for r in leftover_rules:
    if r in list(snakemake_rule_dict.keys()):
        dest = os.path.join(DEST, os.path.relpath(snakemake_rule_dict[r], SNAKEMAKE_RULES_PATH))
        sync_file(snakemake_rule_dict[r], dest, args.dry_run)
    else:
        logger.warn("No such rule '{}' in snakemake_rule_dict".format(r))
