# Copyright (c) 2014 Per Unneberg
# Modelled on bokeh setup script
# --------------------------------------------------
# Imports
# --------------------------------------------------

from __future__ import print_function

# stdlib
import os
from setuptools import setup
from os.path import realpath, dirname, relpath, join

# Extensions
import versioneer

# --------------------------------------------------
# globals and constants
# --------------------------------------------------

ROOT = dirname(realpath(__file__))

# --------------------------------------------------
# classes and functions
# --------------------------------------------------

package_data = []

def package_path(path, filters=()):
    if not os.path.exists(path):
        raise RuntimeError("packaging non-existent path: %s" % path)
    elif os.path.isfile(path):
        package_data.append(relpath(path, 'snakemake_rules'))
    else:
        for path, dirs, files in os.walk(path):
            path = relpath(path, 'snakemake_rules')
            for f in files:
                if not filters or f.endswith(filters):
                    package_data.append(join(path, f))

rule_suffixes = ('.rules', '.rule', '.settings')
                    
package_path(join(ROOT, 'snakemake_rules'), rule_suffixes)
package_path(join(ROOT, 'tests', 'examples', 'Snakefile'))
package_path(join(ROOT, 'tests', 'examples', 'Snakefile_regions'))
package_path(join(ROOT, 'tests', 'examples', 'config.yaml'))
package_path(join(ROOT, 'tests', 'examples', 'config_regions.yaml'))
package_path(join(ROOT, 'tests', 'data'))

scripts = ["scripts/syncrules.py"]

REQUIRES = [
    'snakemake>=3.9.0',
    'pytest',
    'pytest-runner',
]

_version = versioneer.get_version()
_cmdclass = versioneer.get_cmdclass()

setup(
    name="snakemake-rules",
    version=_version,
    cmdclass=_cmdclass,
    author="Per Unneberg",
    author_email="per.unneberg@scilifelab.se",
    description="Snakemake rule library",
    license="MIT",
    url="http://github.com/percyfal/snakemake-rules",
    scripts=scripts,
    packages=[
        'snakemake_rules',
    ],
    package_data={'snakemake_rules': package_data},
    install_requires=REQUIRES,
    setup_requires=['pytest-runner'],
    tests_require=["pytest"],
)
