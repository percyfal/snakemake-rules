Quickstart
==========

Introduction
------------

The purpose of snakemake-rules is to build a library of rules that
can be reused without actually writing them anew. The motivation is
that only parameters, e.g. program options, inputs and outputs, of a
rule change from time to time, but the rule execution is identical.

Installation (optional)
------------------------------

Install using the `conda package manager
<http://conda.pydata.org/docs/>`_:

.. code-block:: shell

   conda install -c percyfal snakemake-rules

Alternatively, install directly from github with pip using the
following command:

.. code-block:: shell
		
   pip3 install -e git+https://github.com/percyfal/snakemake.git@master#egg=snakemake --user


Getting started
---------------

Most importantly, snakemake-rules offers a library of rules that
can be included in a Snakefile and configured via an external
configuration file. Snakemake works in a manner similar to `GNU Make
<https://www.gnu.org/software/make/>`_ in that rules determine how to
generate output files from input files. Output file names are matched
against input file names, whereby wildcards can be used to write
general rules. This feature has been adopted heavily in
snakemake-rules. In fact, most rules are of the form

.. code:: python

   rule somerule:
       input: "{prefix}.inputsuffix"
       output: "{prefix}.outputsuffix"
       run: "command {input} > {output}"

where ``{prefix}`` denotes a wildcard. ``.inputsuffix`` and
``.outputsuffix`` are generally application dependent. For instance,
in the bwa example that follows, the input suffix is generally
``.fastq.gz`` and the output suffix ``.bam``.

If you have installed snakemake-rules, creating a `Snakefile
<https://bitbucket.org/johanneskoester/snakemake/wiki/Documentation#markdown-header-writing-snakefiles>`_
with the following content:

.. code:: python

   # -*- snakemake -*-
   import os
   from snakemake-rules import SNAKEMAKE_RULES_PATH

   include: os.path.join(SNAKEMAKE_RULES__PATH, "bio/ngs/align/bwa.rules")

will add the rules in `bwa.rules`. If you haven't installed
snakemake-rules, you can include the files via urls:

.. code:: python

   # -*- snakemake -*-
   include: "https://raw.githubusercontent.com/percyfal/snakemake-rules/master/snakemake_rules/bio/ngs/align/bwa.rules"

In either case, running

.. code:: shell

	  snakemake -l

should produce

.. code:: shell
	  
   bwa_index
	  bwa index a reference
   samtools_index
	  Run samtools index
   bwa_mem
	  Run bwa mem

To actually run the rules, we need to use the configuration utilities
of snakemake.
