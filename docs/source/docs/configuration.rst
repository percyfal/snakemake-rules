Configuration guide
====================

The configuration guide provides a basic overview of the general
options and settings layout. Snakemake can be configured through a
configuration file that is passed either via the ``--configfile``
command line option, or the Snakefile ``configfile:`` directive. Once
loaded, the configuration settings can be accessed through the global
python object ``config``.

Internally, configuration objects are `python dictionaries`_, where
the keys correspond to configuration options. This has the unfortunate
consequence that it is difficult to provide a documentation API to the
options. This text tries to address this issue, albeit in an
unsufficient manner. As a last resort, for now at least, one simply
has to look at the source code to get an idea of what the options do.
In most cases though, the key names themselves should give an idea of
what behaviour they target.

It is important to keep in mind that no validation of user-supplied
configuration files is done. Consequently, should the user supply a
non-defined configuration key, it will pass unnoticed by Snakemake.
This can be frustrating when debugging; you are sure that you have
changed a configuration value, only to notice later that the
configuration key was misspelled.

Implementation
^^^^^^^^^^^^^^^^

The configuration is constructed as a hierarchy of at most three
levels:

.. code-block:: yaml

   section:
     subsection:
       option:


The ``section`` level corresponds to an application, or a configuration
group of more general nature. The ``subsection`` can either be a new
configuration grouping, or an option to be set. For applications, the
``subsection`` often corresponds to a given rule. Finally, at the
``option`` level, an option is set.

Configuration sections
^^^^^^^^^^^^^^^^^^^^^^

In the toplevel snakemake_rules directory, there is a subdirectory
named ``rules``. The directory contains rules organized by application
directories and ``settings`` files that provide default configuration
values. Every rule directory has its own settings file. There are two
top-level settings file located directly in the ``rules`` directory,
namely ``main.settings`` and ``ngs.settings``.


settings
~~~~~~~~~

The ``settings`` section defines configurations of a general nature. 

.. code-block:: yaml

   settings:
     sampleinfo: sampleinfo.csv
     email: # email
     java:
	java_mem: 8g
	java_tmpdir: /tmp
     runfmt: "{SM}/{SM}_{PU}"
     samplefmt: "{SM}/{SM}"
     threads: 8
     temporary_rules:
       - picard_merge_sam
	
For all settings, see ``rules/main.settings``.

Importantly, many of these settings are *inherited* by the application
rules, so that changing ``threads`` to 4 in ``settings``, will set the
number of threads for all configurations that inherit this option.
However, you can fine-tune the behaviour of the inheriting rules to
override the value in ``settings``; see :ref:`application-settings`.

Here, the most important option is ``sampleinfo``, which provides
information about samples and their naming schemes. For anything but
the smallest workflows, you almost surely need to define this file.
The ``runfmt`` and ``samplefmt`` options describe how the data is
organized. They represent `python miniformat strings`_, where the
entries correspond to columns in the sampleinfo file; hence, in this
case, the column **SM** and **PU** must be present in the sampleinfo
file. So, given the following sampleinfo file

.. code-block:: text

   SM,PU,DT,fastq
   s1,AAABBB11XX,010101,s1_AAABBB11XX_010101_1.fastq.gz
   s1,AAABBB11XX,010101,s1_AAABBB11XX_010101_2.fastq.gz
   s1,AAABBB22XX,020202,s1_AAABBB22XX_020202_1.fastq.gz
   s1,AAABBB22XX,020202,s1_AAABBB22XX_020202_2.fastq.gz

samplefmt will be formatted as ``s1/s1`` and runfmt as
``s1/s1_AAABBB11XX`` or ``s1/s1_AAABBB22XX``, depending on the run. The
formatted strings are used in the workflows as *prefixes* to identify
targets. Rules that operate on the runfmt will be prefixed by
``s1/s1_AAABBB11XX`` or ``s1/s1_AAABBB22XX``, rules that operate on the
sample level (i.e. after merging) will be prefixed by ``s1/s1``.


ngs.settings
~~~~~~~~~~~~~

.. warning::

   The ngs.settings section is slightly disorganized.

``ngs.settings`` affect settings related to ngs analyses:

.. code-block:: yaml
   
   ngs.settings:
     annotation:
	   annot_label: ""
	   transcript_annot_gtf: "",
	   sources: []
     db:
	   dbsnp: ""
       ref: ref.fa
	   transcripts: []
	   build: ""
     fastq_suffix: ".fastq.gz"
     read1_label: "_1"
     read2_label: "_2"
     read1_suffix: ".fastq.gz"
     read2_suffix: ".fastq.gz"
     regions: []
     sequence_capture:
       bait_regions: []
	   target_regions: []


For all settings, see ``rules/ngs.settings``.



samples and ignore_samples
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ``samples`` and ``ignore_samples`` sections are two of the few
top-level configuration keys that are actually set. ``samples`` is a
list of sample names, as defined by the **SM** column in the
``sampleinfo`` file. If no samples are listed, by default all samples
in ``sampleinfo`` will be analyzed. Alternatively, samples can be
ignored through the ``ignore_samples`` key.

.. _application-settings:

Application settings
~~~~~~~~~~~~~~~~~~~~~~~~~

Applications, i.e. bioinformatics software, are grouped in sections by
their application name. Subsections correspond to rules, or
subprograms. For instance, the entire bwa section looks as follows
(with a slight abuse of notation as we here mix yaml with python
objects):

.. code-block::  yaml
   
   bwa:
     cmd: bwa
     ref: config['ngs.settings']['db']['ref']
     index: ""
     index_ext: ['.amb', '.ann', '.bwt', '.pac', '.sa']
     threads: config['settings']['threads']
     mem:
       options:
  

Setting option ``threads`` would then override the value in ``settings``,
providing a means to fine-tune options on a per-application basis.


User-defined configuration
~~~~~~~~~~~~~~~~~~~~~~~~~~

User-defined configurations can be defined in the Snakefile through
the ``config`` object, or through the ``configfile`` directive.
Remember that the configurations must be set **before** any
``include`` statement are issued. The reason is that when a rules file
is included, the default configuration values are compared to the
existing ``config``. If the user has defined custom configurations,
these will take precedence over the default values. If no custom
configuration exists, the default values are applied.


Example
~~~~~~~~

The following example shows two ways to configure the number of
``threads`` and ``options`` for the ``mem`` rule. First, one can
initialize the config object in the Snakefile:

.. code-block:: python

   from os.path import join
   from snakemake_rules import SNAKEMAKE_RULES_PATH
   
   config = {
       'bwa': {
           'mem': {
	       'options': "-M",
	    }
	    'threads': 2,
	}
   }
   include: join(SNAKEMAKE_RULES_PATH, "bwa/bwa_mem.rule")


Alternatively, the configuration can be put in a yaml-file (e.g.
``config.yaml``):

.. code-block:: yaml

   bwa:
     mem:
       options: -M
     threads: 2


and loaded in the Snakefile via the ``configfile`` directive:

.. code-block:: python

   from os.path import join
   from snakemake_rules import SNAKEMAKE_RULES_PATH
   
   configfile: "config.yaml"

   include: join(SNAKEMAKE_RULES_PATH, "bwa/bwa_mem.rule")


.. _python dictionaries: https://docs.python.org/3.5/tutorial/datastructures.html#dictionaries
.. _python miniformat strings: https://docs.python.org/3/library/string.html#formatspec
