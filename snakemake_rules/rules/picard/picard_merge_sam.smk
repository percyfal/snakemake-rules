# -*- snakemake -*-
include: "picard.settings.smk"
include: "picard_sort_sam.smk"

def _picard_merge_sam_input_fn(wildcards):
    """Function to find input files for picard merge.

    Utilise information in config['settings']['sampleinfo'],
    config['settings']['runfmt'] and config['settings']['samplefmt'].
    'runfmt' and 'samplefmt' are python miniformat strings. The
    formatted runfmt string is used to generate names of input bam
    files. The current wildcard prefix corresponds to a formatted
    samplefmt string, all of which are computed and compared to the
    former.
    """
    try:
        suffix = config['picard']['merge_sam']['source_suffix']
        runfmt = config['settings']['runfmt']
        samplefmt = config['settings']['samplefmt']
        samplepfx = {k:k for k in set(samplefmt.format(**s) for s in config['_sampleinfo'])}
        try:
            current_sample = samplepfx[wildcards.prefix]
        except KeyError:
            raise Exception("no such sample prefix {}; sampleinfo file malformatted?".format(wildcards.prefix))
        inputset = list(set(runfmt.format(**x) + suffix for x in config['_sampleinfo'] if samplefmt.format(**x) == current_sample))
        return inputset
    except:
        raise RuleException("""
        picard_merge_sam: Default merge function failed. Make sure the source
        suffix config['picard']['merge_sam']['source_suffix'] is
        correct. Else, implement function for generating source files and
        set configuration
        config['picard']['merge_sam']['inputfun'] to the
        function or a list of input files. Alternatively, use
        workflow.rules.picard_merge_sam.set_input(function) to set the
        input. The function should return a list of bam files to
        be merged""")

config_default = {'picard' :{'merge_sam' : _picard_config_rule_default.copy()}}
config_default['picard']['merge_sam'].update(
    {
        'inputfun' : _picard_merge_sam_input_fn,
        'options' : "CREATE_INDEX=true",
        'runtime' : '24:00:00',
        'source_suffix' : ".bam",
    })


update_config(config_default, config)
config = config_default


rule picard_merge_sam:
    """Picard: merge sam files.

    NB: always outputs bam files!
    """
    params: cmd = config['picard']['cmd'] + MERGE_SAM_FILES,
            options = " ".join([config['picard']['options'],
                                config['picard']['merge_sam']['options']]),
            runtime = config['picard']['merge_sam']['runtime']
    input: config['picard']['merge_sam']['inputfun']
    output: merge="{prefix}.merge.bam"
    threads: config['picard']['merge_sam']['threads']
    run:
      if (len(input) > 1):
          inputstr = " ".join(["INPUT={}".format(x) for x in input])
          shell("{cmd} {ips} OUTPUT={out} {opt}".format(cmd=params.cmd, ips=inputstr, out=output.merge, opt=params.options))
      else:
          if os.path.exists(output.merge):
              os.unlink(output.merge)
          shell("{picard} {cmd} I={infile} O={out} {opt} {sortopt}".format(picard=config['picard']['cmd'],
                                                                           cmd=SORT_SAM,
                                                                           infile=input[0], out=output.merge,
                                                                           opt=config['picard']['options'],
                                                                           sortopt=config['picard']['sort_sam']['options']))
