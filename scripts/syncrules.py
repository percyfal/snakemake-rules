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
    if not os.path.exists(v):
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
rules = {x:os.path.join(DEST, os.path.relpath(x, SNAKEMAKE_RULES_PATH)) for x in workflow.included if x.startswith(SNAKEMAKE_RULES_PATH)}

# Copy rules to outdir
for k, v in rules.items():
    sync_file(k, v, args.dry_run)
