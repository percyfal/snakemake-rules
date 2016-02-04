Configuration
=============

Rule organization
~~~~~~~~~~~~~~~~~~~

Rules are organized by application in *rule application files* that
are grouped by a multi-level directory structure that reflects the
application type. For instance, the file ``bwa.rules`` for the
alignment application `bwa <http://bio-bwa.sourceforge.net/>`_ is
located in the directory ``bio/ngs/align``. Every application has a
rule application file that ends in the suffix ``.rules``.

However, the rule application file does not contain any actual rules.
Its main purpose is to load and initialize global configuration
variables, and to define default configuration values applicable to
*all* rules for the given application. The actual application rules
are stored one rule per file in a subdirectory named with an
underscore and the application name. Thus, the subdirectory ``_bwa``
stores rules for bwa in files named after the rule ending in the
suffix ``.rule``.

It is important to remember that the ``.rule`` files should never be
directly included in Snakefiles. Which rules are included in defined
by a configuration parameter in the application rule file called
``rules``. For instance, an excerpt from ``bwa.rules`` looks as
follows:

.. code-block:: python

   include: '../settings.rules'

   DEFAULT_RULES = ['bwa_index', 'bwa_mem']

   config_default = {
      'bio.ngs.align.bwa' : {
          ...
	  'rules' : DEFAULT_RULES,
	  },
   }

   update_config(config_default, config)
   config = config_default

   for rule in config['bio.ngs.align.bwa']['rules']:
       include: os.path.join("_bwa", rule + ".rule")

The rules ``bwa_index`` and ``bwa_mem`` are loaded (included) by
default, but by modifying the ``rules`` parameter, one can specify what
rules to include. The following sections will show more specifically
how this is done.

The default configuration
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For any application, the application rule file contains a *default
configuration* that ensures that all rules have sensible defaults set,
regardless whether the user decides to modify them or not. The basic
configuration structure looks like

.. code-block:: text

   namespace
       section/parameter
           parameter

The *namespace* is an identifier for the rules file, and should be
named ``path.to.rules``, where ``path`` and ``to`` are directory names
relative to the rules root path. The *section/parameter* is either a
parameter related to the program, or a subprogram which in turn can
have *parameters* assigned to it. The section levels are mapped to the
corresponding levels in the snakemake global *config* dictionary
object.

As an example, consider the default configuration for
bwa once again:

.. code-block:: python

   config_default = {
       'bio.ngs.align.bwa' : {
           'cmd' : "bwa",
	   'ref' : config['bio.ngs.settings']['db']['ref'],
	   'index' : "",
	   'index_ext' : ['.amb', '.ann', '.bwt', '.pac', '.sa'],
           'threads' : config['settings']['threads'],
           'rules' : DEFAULT_RULES,
       },
   }
   update_config(config_default, config)
   config = config_default

   
The namespace is ``bio.ngs.align.bwa``, reflecting the fact that the
application rules file is located in the folder ``bio/ngs/align`` and
is named ``bwa.rules``. ``config`` is the global snakemake
configuration object. The ``config_default`` dictionary is mapped to
the ``config`` object by use of the function
:meth:`snakemake.utils.update_config`. First, ``config_default`` is
updated with the current state of ``config``. Should ``config`` define
any of the parameters in ``config_default``, the latter are
overwritten; if not, the default values are used. The updated
``config_default`` dictionary is then set to ``config``. After the
configuration has been updated, the ``threads`` parameter above would
be accessible via ``config['bio.ngs.align.bwa']['threads']``.

Incidentally, this example shows another key idea of the
configuration, namely that some options inherit from rules files
higher up in the file hierarchy. The rules file
``bio/ngs/settings.rules`` contains a generic configuration that is
common to all ngs rules. This implementation makes it possible to
override settings for specific programs, like for instance in the case
of the ``threads`` parameter.

The application configuration is applicable to all bwa rules.
However, each bwa rule can have additional configuration defaults
that are defined in the rule file itself. The top section of
``bwa_mem.rule`` reads as follows:

.. code-block:: python

   config_default = {
       'bio.ngs.align.bwa': {
           'mem' : {
               'options': "",
	   },
       },
   }



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
       'bio.ngs.align.bwa': {
           'mem': {
	       'options': "-M",
	    }
	    'threads': 2,
	}
   }
   include: join(SNAKEMAKE_RULES_PATH, "bio/ngs/align/bwa.rules")


Alternatively, the configuration can be put in a yaml-file (e.g.
``config.yaml``):

.. code-block:: yaml

   bio.ngs.align.bwa:
     mem:
       options: -M
     threads: 2


and loaded in the Snakefile via the ``configfile`` directive:

.. code-block:: python

   from os.path import join
   from snakemake_rules import SNAKEMAKE_RULES_PATH
   
   configfile: "config.yaml"

   include: join(SNAKEMAKE_RULES_PATH, "bio/ngs/align/bwa.rules")
