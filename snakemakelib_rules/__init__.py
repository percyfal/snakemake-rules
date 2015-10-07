from ._version import get_versions
__version__ = get_versions()['version']
del get_versions

# Set rule path
import os
SNAKEMAKELIB_RULES_PATH = os.path.dirname(__file__)
