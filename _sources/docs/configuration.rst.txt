Configuration
=============

Rule organization
~~~~~~~~~~~~~~~~~~~

Rules are organized by application directories. Each directory
contains a *settings* file, that initializes global configuration
variables, and to define default configuration values applicable to
*all* rules for the given application. The actual application rules
are stored one rule per file with suffix ``.rule``. For instance, the
file ``bwa_mem.rule`` for the alignment application `bwa
<http://bio-bwa.sourceforge.net/>`_ contains a rule for the command
``bwa mem`` and is located in the ``bwa`` directory.


The default configuration
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For any rule, there is a *default configuration* that ensures that
defaults are set, regardless whether the user decides to modify them
or not. The basic configuration structure looks like

.. code-block:: text

   application
       section/parameter
           parameter

The *application* is an identifier for the application directory. The
*section/parameter* is either a parameter related to the program, or a
subprogram which in turn can have *parameters* assigned to it. The
section levels are mapped to the corresponding levels in the snakemake
global *config* dictionary object.

As an example, consider the default configuration for
bwa once again:

.. code-block:: python

   config_default = {
       'bwa' : {
           'cmd' : "bwa",
	   'ref' : config['ngs.settings']['db']['ref'],
	   'index' : "",
	   'index_ext' : ['.amb', '.ann', '.bwt', '.pac', '.sa'],
           'threads' : config['settings']['threads'],
       },
   }
   update_config(config_default, config)
   config = config_default

   
The application is ``bwa``, reflecting the fact that the rule files
 are located in the folder ``bwa``. ``config`` is the global snakemake
 configuration object. The ``config_default`` dictionary is mapped to
 the ``config`` object by use of the function
 :meth:`snakemake.utils.update_config`. First, ``config_default`` is
 updated with the current state of ``config``. Should ``config``
 define any of the parameters in ``config_default``, the latter are
 overwritten; if not, the default values are used. The updated
 ``config_default`` dictionary is then set to ``config``. After the
 configuration has been updated, the ``threads`` parameter above would
 be accessible via ``config['bwa']['threads']``.

Incidentally, this example shows another key idea of the
configuration, namely that some options inherit from rules files
higher up in the file hierarchy. The settings file ``main.settings``
contains a generic configuration that is loaded by all application
settings. This implementation makes it possible to override settings
for specific programs, like for instance in the case of the
``threads`` parameter.

The application configuration is applicable to all bwa rules.
However, each bwa rule can have additional configuration defaults
that are defined in the rule file itself. The top section of
``bwa_mem.rule`` reads as follows:

.. code-block:: python

   config_default = {
       'bwa': {
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
