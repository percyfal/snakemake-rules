Examples
==========

Bwa alignment
---------------

.. note:: This example requires you have installed `bwa
          <http://bio-bwa.sourceforge.net/>`_ and `samtools
          <http://www.htslib.org/>`_

First, create a `Snakefile
<https://bitbucket.org/johanneskoester/snakemake/wiki/Documentation#markdown-header-writing-snakefiles>`_
with the following content:

.. code:: python

   # -*- snakemake -*-
   import os
   from snakemakelib-rules import SNAKEMAKELIB_PATH
   
   config = {
       'bio.ngs.align.bwa' : {
	   'index' : os.path.join(SNAKEMAKELIB_PATH, "data/genomes/Hsapiens/hg19/bwa/chr11.fa"),
       },
   }

   workdir: os.path.join(SNAKEMAKELIB_PATH, "data/projects/J.Doe_00_01")
   include: os.path.join(SNAKEMAKELIB_PATH, "rules/bio/ngs/align/bwa.rules")

Briefly, in this file we:

1. set a **configuration dictionary** that indicates where the bwa
   index files are located
2. set a **working directory**
3. **include** rules for bwa

In this minimal example, we make use of the internal snakemakelib
variable ``SNAKEMAKELIB_PATH`` to locate the rules, but could as well
have included the rules by supplying the full path.

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

Now, in the test data set relative to the working directory defined
above there are sequence input files organized under samples and
sequencing runs. For instance, in subfolder
``P001_101/120924_AC003CCCXX`` we have the files

.. code:: shell

   1_120924_AC003CCCXX_P001_101_1.fastq.gz
   1_120924_AC003CCCXX_P001_101_2.fastq.gz

Since the bwa rule has the output suffix ``.bam`` and input suffix
``.fastq.gz``, we can align these input files by issuing [#f1]_

.. code:: shell

   snakemake -F P001_101/120924_AC003CCCXX/1_120924_AC003CCCXX_P001_101.bam

In addition to performing the alignment, this command will generate
bwa indices on the fly. The flag ``-F`` tells snakemake to rerun the
rules, even if outputs are present.

