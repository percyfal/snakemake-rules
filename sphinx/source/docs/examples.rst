==========
 Examples
==========

snakemake-rules is shipped with some simple tests, including test data
files and Snakefiles with auxiliary rules to run simple pipelines. The
following examples are based on the example Snakefiles, located in the
``/tests/examples`` subdirectory of the ``snakemake_rules``
installation directory . In the examples, we will use yaml
configuration files that are loaded via the command line.


Bwa alignment
=============

.. note:: To actually run this example requires `bwa
          <http://bio-bwa.sourceforge.net/>`_ and `samtools
          <http://www.htslib.org/>`_.

First, cd to the test directory. It contains a `Snakefile
<https://bitbucket.org/johanneskoester/snakemake/wiki/Documentation#markdown-header-writing-snakefiles>`_
called ``Snakefile_regions`` with the following content:

.. code:: python

   # -*- snakemake -*-
   from os.path import join
   from snakemake_rules import SNAKEMAKE_RULES_PATH

   include: join(SNAKEMAKE_RULES_PATH, "bwa/bwa_index.rule")
   include: join(SNAKEMAKE_RULES_PATH, "bwa/bwa_mem.rule")
   include: join(SNAKEMAKE_RULES_PATH, "samtools/samtools_sort.rule")
   include: join(SNAKEMAKE_RULES_PATH, "samtools/samtools_index.rule")
   
Briefly, in this file we **include** rules for bwa and samtools. In
this minimal example, we make use of the internal snakemake variable
``SNAKEMAKE_RULES_PATH`` to locate the rules, but could as well have
included the rules by supplying the full path. By setting the working
directory via the command line, we can place the Snakefile anywhere we
like; the commands will still be executed in the working directory.

In the same directory there is also the configuration file
``config_regions.yaml`` with the following content:

.. code-block:: yaml

   bwa:
     index: ref.fa

Here, we indicate that the bwa index file ``ref.fa`` is located in the
working directory. 

Now, to see which rules are included, you can type:

.. code:: shell

   $ snakemake -l

which should generate the following output:
   
.. code:: shell

   bwa_index
       bwa index a reference
   bwa_mem
       Run bwa mem
   samtools_sort
       Run samtools sort
   samtools_index
       Run samtools index
   all
    
The bwa rule has the output suffix ``.bam`` and input suffix
``.fastq.gz``. Therefore, we can align the input files by issuing

.. code:: shell

   $ snakemake -d ../data --configfile config_regions.yaml -s Snakefile_regions

Prior to running a command, it is advisable to use the flags ``-n``
(dry run) and ``-p`` (print commands) to actually see what will
happen. Optionally, one can add ``-F`` which will force snakemake to
rerun rules, even if output files should exist. Then,

.. code:: shell

   $ snakemake -s Snakefile_regions -n -p -F s1.bam --configfile config_regions.yaml -d ../data

generates the output

.. code:: shell

   rule bwa_index:
       input: ref.fa
       output: ref.fa.amb, ref.fa.ann, ref.fa.bwt, ref.fa.pac, ref.fa.sa
       wildcards: prefix=ref, ext=.fa

   bwa index ref.fa
   rule bwa_mem:
       input: ref.fa.amb, ref.fa.ann, ref.fa.bwt, ref.fa.pac, ref.fa.sa, s1_2.fastq.gz, s1_1.fastq.gz
       output: s1.bam
       log: s1.log
       wildcards: prefix=s1

   bwa mem -t 1  ref.fa s1_1.fastq.gz s1_2.fastq.gz 2> s1.log |  samtools view -Sb - > s1.bam
   Job counts:
	   count	jobs
	   1	bwa_index
	   1	bwa_mem
	   2

Consequently, ``bwa index`` will first be run on ``ref.fa``, followed
by ``bwa mem`` on the input sequence files, and so on.


Variant calling
===============

A slightly more complicated example is given in the ``Snakefile``.
However, the only major difference to the previous example is that
more rules have been included, and there is an  ``all`` rule:

.. code-block:: python

   rule all:
       input: "s1.sort.rg.dup.bcftools.vcf.gz"

By combining suffixes in the right order and defining a desired output
file, we generate a pipeline on the fly. The workflow can be
visualized with the command

.. code-block:: shell

   $ snakemake -n -p --rulegraph

with the following result

.. graphviz::

   digraph snakemake_dag {
       graph[bgcolor=white, margin=0];
       node[shape=box, style=rounded, fontname=sans,                 fontsize=10, penwidth=2];
       edge[penwidth=2, color=grey];
	   0[label = "bwa_mem", color = "0.40 0.6 0.85", style="rounded"];
	   1[label = "samtools_sort", color = "0.13 0.6 0.85", style="rounded"];
	   2[label = "all", color = "0.00 0.6 0.85", style="rounded"];
	   3[label = "bam_fofn", color = "0.20 0.6 0.85", style="rounded"];
	   4[label = "picard_add_or_replace_read_groups", color = "0.53 0.6 0.85", style="rounded"];
	   5[label = "bwa_index", color = "0.07 0.6 0.85", style="rounded"];
	   6[label = "bcftools_call", color = "0.33 0.6 0.85", style="rounded"];
	   7[label = "picard_build_bam_index", color = "0.27 0.6 0.85", style="rounded"];
	   8[label = "picard_mark_duplicates", color = "0.60 0.6 0.85", style="rounded"];
	   5 -> 0
	   0 -> 1
	   6 -> 2
	   8 -> 3
	   1 -> 4
	   3 -> 6
	   4 -> 7
	   7 -> 8
	   4 -> 8
   }                 

Merging files
=============

Some rules require some additional tinkering; ``picard_merge_sam`` is
one such rule. It falls into the class of rules that depend on input
files whose number and names need to be generated by the rule itself.
To see what this means in practice, let's take a closer look at what
the rule does.

The `picard <http://broadinstitute.github.io/picard/>`_ command
`MergeSamFiles
<http://broadinstitute.github.io/picard/command-line-overview.html#MergeSamFiles>`_
merges several SAM/BAM files into one file. The input files should be
provided as a list to snakemakes ``input`` directive. Instead of
hard-coding names, Snakemake lets ``input`` take a `function as
argument
<https://bitbucket.org/snakemake/snakemake/wiki/Documentation#markdown-header-functions-as-input-files>`_,
a function that **must** take one parameter ``wildcards``. In the case
of ``picard_merge_sam``, the function should return the list of files
to be merged. Note that the function cannot simply look for existing
output files (unless of course they are the starting files for the
analysis); the rule needs to generate file names that are produced by
some upstream rule.

The ``input`` and ``output`` section of ``picard_merge_sam`` is shown
below:

.. code-block:: python

   input: config['picard']['merge_sam']['inputfun']
   output: merge="{prefix}." + config['picard']['merge_sam']['label'] + ".bam"


The configuration entry
``config['picard']['merge_sam']['inputfun']`` should be set
to the function in question. There is one wildcard ``prefix`` which in
the function is accessible through ``wildcards.prefix``.

Suppose for instance we have a merge target ``s.merge.bam`` that
takes as input files ``s1.bam`` and ``s2.bam``. Then the
following python code defines a function that generates the correct
file names and sets the relevant configuration section:

.. code-block:: python

   def merge_inputs(wildcards):
       return [wildcards.prefix + "1.bam", wildcards.prefix + "2.bam"] 


   config = {'picard': {'merge_sam': {'inputfun': merge_inputs}}}
