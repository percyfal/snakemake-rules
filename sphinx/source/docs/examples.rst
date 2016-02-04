Examples
==========

snakemake_rules is shipped with some simple tests, including test data
files and Snakefiles with auxiliary rules to run simple pipelines. The
following examples are based on the example Snakefiles, located in the
directory ``snakemake_rules/tests``. In the examples, we will use yaml
configuration files that are loaded via the ``configfile`` directive.


Bwa alignment
---------------

.. note:: To actually run this example requires `bwa
          <http://bio-bwa.sourceforge.net/>`_ and `samtools
          <http://www.htslib.org/>`_

First, create a `Snakefile
<https://bitbucket.org/johanneskoester/snakemake/wiki/Documentation#markdown-header-writing-snakefiles>`_
with the following content:

.. code:: python

   # -*- snakemake -*-
   from os.path import join
   from snakemake_rules import SNAKEMAKE_RULES_PATH

   workdir: join(SNAKEMAKE_RULES_PATH, "tests")

   configfile: "config_bwa.yaml"
   
   include: join(SNAKEMAKE_RULES_PATH, "bio/ngs/align/bwa.rules")

   
Briefly, in this file we set a **working directory**, load the
**configuration file** and **include** rules for bwa. In this minimal
example, we make use of the internal snakemake variable
``SNAKEMAKE_RULES_PATH`` to locate the rules, but could as well have
included the rules by supplying the full path. By setting the working
directory here, we can place the Snakefile anywhere we like; the
commands will still be executed in the working directory.

Create the configuration file ``config_bwa.yaml`` with the following
content:

.. code-block:: yaml

   bio.ngs.bwa.align:
     index: data/chr11.fa

Here, we indicate that the bwa index file ``chr11.fa`` is located in a
directory ``data`` relative to the working directory. The ``data``
directory also contains a pair of sequence input files,
``test_1.fastq.gz`` and ``test_2.fastq.gz``.

Now, to see which rules are included, you can type:

.. code:: shell

   $ snakemake -l

which should generate the following output:
   
.. code:: shell

   bwa_index
	bwa index a reference
   samtools_index
	Run samtools index
   bwa_mem
	Run bwa mem

So, by including ``bwa.rules``, we have actually defined three
`snakemake rules
<https://bitbucket.org/johanneskoester/snakemake/wiki/Documentation#markdown-header-rules>`_.
The bwa rule has the output suffix ``.bam`` and input suffix
``.fastq.gz``. Therefore, we can align the input files by issuing

.. code:: shell

   $ snakemake data/test.bam

Prior to running a command, it is advisable to use the flags ``-n``
(dry run) and ``-p`` (print commands) to actually see what will
happen. Optionally, one can add ``-F`` which will force snakemake to
rerun rules, even if output files should exist. Then,

.. code:: shell

   $ snakemake -n -p -F data/test.bam

generates the output

.. code:: shell

   rule bwa_index:
	input: data/chr11.fa
	output: data/chr11.fa.amb, data/chr11.fa.ann, data/chr11.fa.bwt, data/chr11.fa.pac, data/chr11.fa.sa
   bwa index data/chr11.fa
   rule bwa_mem:
	input: data/chr11.fa.amb, data/chr11.fa.ann, data/chr11.fa.bwt, data/chr11.fa.pac, data/chr11.fa.sa, data/test_1.fastq.gz, data/test_2.fastq.gz
	output: data/test.bam
	log: data/test.log
   bwa mem -t 1  data/chr11.fa data/test_1.fastq.gz data/test_2.fastq.gz | samtools view -Sb - > data/test.bam
   Job counts:
           count    jobs
	   1	    bwa_index
	   1	    bwa_mem
	   2

Consequently, ``bwa index`` will first be run on ``data/chr11.fa``,
followed by ``bwa mem`` on the input sequence files.


Variant calling
----------------

A slightly more complicated example is given in the Snakefile in
tests. However, the only major difference to the previous example is
that more applications (rules files) have been included, and a rule
``all`` has been added:

.. code-block:: python

   rule all:
       input: "data/test.sort.rg.dup.realign.recal.annotated.vcf"

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
	0[label = "bwa_mem", color = "0.58 0.6 0.85", style="rounded"];
	1[label = "bwa_index", color = "0.19 0.6 0.85", style="rounded"];
	2[label = "picard_mark_duplicates", color = "0.59 0.6 0.85", style="rounded"];
	3[label = "gatk_unified_genotyper", color = "0.41 0.6 0.85", style="rounded"];
	4[label = "gatk_realigner_target_creator", color = "0.44 0.6 0.85", style="rounded"];
	5[label = "gatk_print_reads", color = "0.06 0.6 0.85", style="rounded"];
	6[label = "picard_add_or_replace_read_groups", color = "0.43 0.6 0.85", style="rounded"];
	7[label = "snpeff_annotate_variants", color = "0.24 0.6 0.85", style="rounded"];
	8[label = "all", color = "0.09 0.6 0.85", style="rounded"];
	9[label = "gatk_indel_realigner", color = "0.28 0.6 0.85", style="rounded"];
	10[label = "picard_build_bam_index", color = "0.30 0.6 0.85", style="rounded"];
	11[label = "gatk_base_recalibrator", color = "0.52 0.6 0.85", style="rounded"];
	12[label = "samtools_sort", color = "0.53 0.6 0.85", style="rounded"];
	1 -> 0
	6 -> 2
	10 -> 2
	5 -> 3
	2 -> 4
	10 -> 4
	9 -> 5
	11 -> 5
	12 -> 6
	3 -> 7
	7 -> 8
	2 -> 9
	4 -> 9
	2 -> 10
	6 -> 10
	9 -> 11
	0 -> 12
   }  
       

Merging files
--------------

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

   input: config['bio.ngs.qc.picard']['merge_sam']['inputfun']
   output: merge="{prefix}." + config['bio.ngs.qc.picard']['merge_sam']['label'] + ".bam"


The configuration entry
``config['bio.ngs.qc.picard']['merge_sam']['inputfun']`` should be set
to the function in question. There is one wildcard ``prefix`` which in
the function is accessible through ``wildcards.prefix``.

Suppose for instance we have a merge target ``test.merge.bam`` that
takes as input files ``test.run1.bam`` and ``test.run2.bam``. Then the following python code defines a function that generates the correct file names and sets the relevant configuration section:

.. code-block:: python

   def merge_inputs(wildcards):
       return [wildcards.prefix + ".run1.bam", wildcards.prefix + ".run2.bam"] 


   config['bio.ngs.qc.picard']['merge_sam']['inputfun'] = merge_inputs
