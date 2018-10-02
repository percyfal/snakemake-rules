# -*- snakemake -*-
include: 'mapdamage2.settings.smk'

m = re.search("(?P<opt>-d|--folder)\s+(?P<folder>.*)", config["mapdamage2"]["options"])
if m:
    logger.info("--folder given in option string; removing. Set folder via the wildcard starting with md")
    folder = m.group("folder")
    opt = m.group("opt")
    config["mapdamage2"]["options"] = re.sub(opt, "", config["mapdamage2"]["options"])
    config["mapdamage2"]["options"] = re.sub(folder, "", config["mapdamage2"]["options"])

m = re.search("(?P<nostats>--no-stats)", config["mapdamage2"]["options"])
if m:
    pass
else:
    logger.info("adding --no-stats to option string.")
    config["mapdamage2"]["options"] = config["mapdamage2"]["options"] + " --no-stats"

rescalebam = None
m = re.search("(?P<rescalebam>--rescale-out\s+[\w\.]+)", config["mapdamage2"]["options"])
if m:
    logger.warning("--rescale-out given in option string; replacing with prefix based name")
    rescalebam = m.group("rescalebam")
    config["mapdamage2"]["options"] = re.sub(rescalebam, "", config["mapdamage2"]["options"])

rescale = False
m = re.search("(--rescale\s+)", config["mapdamage2"]["options"])
if rescalebam and not m:
    config["mapdamage2"]["options"] += " --rescale "
    rescale = True
elif m:
    rescale = True

rule mapdamage2:
    """mapdamage2

    Run mapdamage2 without stats.
    """
    version: "0.2"
    params: options = config["mapdamage2"]["options"],
            runtime = config["mapdamage2"]["runtime"]
    wildcard_constraints: folder = "(mapdamage2|md[^ ]*)"
    input: bam = "{prefix}.bam", ref = config["mapdamage2"]["ref"]
    output: log = "{prefix}_{folder}/Runtime_log.txt",
            fragmisincorporation = "{prefix}_{folder}/Fragmisincorporation_plot.pdf",
            misincorporation = "{prefix}_{folder}/misincorporation.txt",
            CtoT = "{prefix}_{folder}/5pCtoT_freq.txt",
            GtoA = "{prefix}_{folder}/3pGtoA_freq.txt",
            dnacomp = "{prefix}_{folder}/dnacomp.txt",
            lgdist = "{prefix}_{folder}/lgdistribution.txt",
            #rescale = "{prefix}.rescaled.bam" if rescale else []
    threads: config["mapdamage2"]["threads"]
    run:
        cmd = "mapDamage " + params.options + " -i " + input.bam + " -r " + input.ref + " -d " + "{}_{}".format(wildcards.prefix, wildcards.folder)
        if rescale:
            cmd = cmd + " --rescale-out " + output.rescale
        shell(cmd)
