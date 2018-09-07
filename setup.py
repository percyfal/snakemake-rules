# Copyright (c) 2014 Per Unneberg
# Modelled on bokeh setup script
# --------------------------------------------------
# Imports
# --------------------------------------------------

from __future__ import print_function

# stdlib
import os
import sys
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


rule_suffixes = ('.smk')

package_path(join(ROOT, 'snakemake_rules'), rule_suffixes)
package_path(join(ROOT, 'tests', 'examples', 'Snakefile'))
package_path(join(ROOT, 'tests', 'examples', 'Snakefile_regions'))
package_path(join(ROOT, 'tests', 'examples', 'config.yaml'))
package_path(join(ROOT, 'tests', 'examples', 'config_regions.yaml'))
package_path(join(ROOT, 'tests', 'data'))

scripts = ["scripts/syncrules.py"]

setup_requires = []
if {'pytest', 'test'}.intersection(sys.argv):
    setup_requires.append('pytest-runner')

install_requires = [
    'snakemake>=5.0.0',
]

test_requires = [
    'pytest',
    'pytest-runner',
    'psutil',
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
        'snakemake_rules.core',
    ],
    package_data={'snakemake_rules': package_data},
    install_requires=install_requires,
    setup_requires=setup_requires,
    tests_require=test_requires,
)
