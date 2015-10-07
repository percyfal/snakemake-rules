.. snakemakelib-rules documentation master file, created by
   sphinx-quickstart on Tue Oct  6 15:00:15 2015.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

snakemakelib-rules - a library of snakemake rules
==================================================

.. _about:


`Snakemake <https://bitbucket.org/johanneskoester/snakemake/wiki/Home>`__
library for various applications, with a focus on bioinformatics and
next-generation sequencing.

The Snakemake rules contain general recipies for commonly used
applications and bioinformatics programs. The use cases reflect the
needs I've had and do by no means have a comprehensive coverage.
Nevertheless, many commands are so commonly used that the recipes may
be of general interest.

.. warning:: Use the rules at your own risk, and make sure you
             understand them before running any commands. I take no
             responsibility if you'd happen to run a ``snakemake
             clean`` in an inappropriate location, removing precious
             data in the process.


Features
^^^^^^^^

1. **Rule library**. snakemakelib-rules is just a library of snakemake
   rules. At the very least, if rules need to be tweaked, the rule
   library serves as a cut-and-paste resource of template rules
2. **Atomic rules**. As much as is possible, every rule lives in an
   individual file. This makes it easy to fine-tune what rules to
   include via application-specific rules file.
3. **Default configuration**. All rules contain configuration defaults
   that populate the snakemake configuration variable `config`. The
   configuration defaults can be overridden with local configuration
   files.
4. **Installation free**. Most rules can be included in snakefiles
   without even using functionality in snakemakelib. You can use the
   github url to include rules:

   .. code-block:: python

      include: "https://raw.githubusercontent.com/percyfal/snakemakelib-rules/master/snakemakelib_rules/bio/ngs/align/bwa.rules"



Contents
---------

.. toctree::
   :maxdepth: 2

   docs/quickstart
   docs/configuration
   docs/examples
   docs/release_notes
