# -*- snakemake -*-
include: "picard.settings.smk"

def _add_rg(wildcards):
    """Return read group string"""
    current_sample = wildcards.prefix
    rg = {'ID': current_sample,
          'SM': current_sample,
          'LIB': 'lib',
          'PU': 'pu',
          'PL': 'Illumina'
    }
    try:
        runfmt = config['settings']['runfmt']
        rginfo = None
        for s in config['_sampleinfo']:
            if current_sample.startswith(runfmt.format(**s)):
                rginfo = s
                continue
        if rginfo is None:
            raise Exception("no such run prefix {}; sampleinfo file malformatted?".format(current_sample))
        rg.update(rginfo)
        rg['rgid'] = "{SM}_{PU}".format(SM=rg['SM'], PU=rg['PU'])
        rg_str = "RGID={ID} RGSM={SM} RGLB={LIB} RGPL={PL} RGPU={PU}".format(**rg)
        logger.info("setting read group information to " + rg_str)
        return rg_str
    except:
        raise RuleException("""Implement function for adding read group information and set
        configuration config['picard']['add_or_replace_read_groups']['rgfun']
        to the function. The function should take a file name as its input parameter.""")


config_default = {'picard' :{'add_or_replace_read_groups' : _picard_config_rule_default.copy()}}
config_default['picard']['add_or_replace_read_groups'].update(
    {
        'options' : "sort_order=coordinate create_index=true",
        'rgfun' : _add_rg,
    })

update_config(config_default, config)
config = config_default

rule picard_add_or_replace_read_groups:
    """Picard: add or replace read groups. Currently tailored for Illumina read groups."""
    params: cmd = config['picard']['cmd'] + ADD_OR_REPLACE_READ_GROUPS,
            options = config['picard']['options'],
            custom_options = config['picard']['add_or_replace_read_groups']['options'],
            runtime = config['picard']['add_or_replace_read_groups']['runtime']
    input: bam = "{prefix}.bam"
    output: bam = "{prefix}.rg.bam"
    threads: config['picard']['add_or_replace_read_groups']['threads']
    run:
        shell("{params.cmd} INPUT={input.bam} OUTPUT={output.bam} {params.options} {params.custom_options} " + config['picard']['add_or_replace_read_groups']['rgfun'](wildcards))
