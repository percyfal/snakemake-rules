# -*- snakemake -*-
import io
from snakemake.io import Wildcards, Namedlist

include: "utils.settings.smk"

rule dag:
    """Print dag

    NB: currently this rule may have to be run manually. When the rule
    is included as a dependecy it fails throwing an
    MissingInputException

    """
    version: "0.1"
    output: dag = "{prefix}_dag.dot"
    run:
        # see http://stackoverflow.com/questions/5136611/capture-stdout-from-a-script-in-python
        backup = sys.stdout
        sys.stdout = io.StringIO()
        workflow.execute(targets=[os.path.basename(wildcards.prefix)],
                         dryrun=False, printdag=True, updated_files=[])
        out = sys.stdout.getvalue()
        sys.stdout.close()
        sys.stdout = backup
        with open (output.dag, "w") as fh:
            fh.write(out)
