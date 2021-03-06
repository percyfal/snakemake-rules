# -*- snakemake -*-
import yaml

rule conf_yaml:
    """Print global configuration to yaml output file"""
    version: "0.1"
    output: yaml = os.path.join("{path}", "config_global.yaml")
    run:
        with open (output.yaml,"w") as fh:
            fh.write(yaml.dump(config, default_flow_style=False))


localrules: conf_yaml
