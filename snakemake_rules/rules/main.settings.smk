# -*- snakemake -*- 
"""Main settings for snakemake-rules. Defines top-level configuration
settings for all rules.
"""
import sys
import os
from os.path import join, basename, dirname
import re
import shutil
import glob
from snakemake_rules.core.utils import get_samples, python2_path, available_mem, cd
from snakemake.utils import update_config
from snakemake.exceptions import RuleException
from snakemake.logging import logger
try:
    from snakemake_rules._version import get_versions
    from snakemake_rules import SNAKEMAKE_RULES_PATH
except:
    get_versions = lambda: {'full-revisionid': 'unknown', 'version': 'unknown', 'dirty': "unknown", 'error': None}

# Check for old configuration key versions
if any(x.startswith("bio.") for x in list(config.keys())):
    raise Exception("bio. namespace is deprecated; please update your configuration keys")

# Get the number of threads passed via command line
_threads = 8
m = re.search("(-j|--jobs|--cores)\s+(?P<threads>\d+)", " ".join(sys.argv))
if m:
    _threads = int(m.group("threads"))

# Set the available memory to memory per core
_mem = available_mem(1, "8g")

config_default = {
    'settings' : {
        'email' : os.getenv("USER_EMAIL", None),
        'java' : {
            'java_mem' : _mem,
            'java_tmpdir' : "/tmp",
        },
        'threads' : _threads,
        'temporary_rules': [],
        'protected_rules': [],
        'sampleinfo': None,
        'conda' : {
            'python2' : 'py2.7',
        }, 
    },
    'samples': {'include': [], 'exclude': []},
}

update_config(config_default, config)
config = config_default

get_samples(config, logger)
python2_path(config, logger)
    
config['_snakemake_rules_version'] = get_versions()
