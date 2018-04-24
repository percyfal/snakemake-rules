# -*- snakemake -*-
include: "bamtools.settings"
include: "bamtools_filter.rule"

import json

rule bamtools_create_filter_script:
    """Create bamtools filter script"""
    params: script = "\n".join('{{\n"reference":"{r}",{s}\n}},'.format(r=r, s=json.dumps(config['bamtools']['filter']['options'], indent=1).lstrip("{").rstrip("}")) for r in config['bamtools']['filter']['regions']).rstrip(",")
    output: script = temporary("{prefix}.script")
    shell:
        """echo '{{"filters":[\n{params.script}\n    ]\n}}\n' > {output.script}"""

localrules: bamtools_create_filter_script        
