# -*- snakemake -*-
# Simplest greedy algorithm: see https://en.wikipedia.org/wiki/Partition_problem#The_greedy_algorithm
# 1. choose the p largest regions to initial partitions
# 2. for each remaining element, add it to the smallest set
include: "pybedtools.settings.smk"

config_default = {
    'pybedtools' : {
        'make_bed_windows' : {
            'window_size' : 100000,
            'partitions' : 1000,
        },
    },
}

update_config(config_default, config)
config = config_default


rule pybedtools_make_bed_windows:
    """Generate bed targets based on window size.
    
    Partition input bed file into windows. Intersect with original bam
    file.

    """
    params: runtime = "01:00:00"
    input: targets = "{prefix}.bed", chrom_sizes = "{prefix}.chrom.sizes"
    output: targets = ["{{prefix}}.{partition}.bed".format(partition="window-{}".format(p+1)) for p in range(config['pybedtools']['make_bed_windows']['partitions'])]
    threads: 1
    run:
        logger.info("loading targets {}".format(input.targets))
        targets = pybedtools.BedTool(input.targets)
        w = config['pybedtools']['make_bed_windows']['window_size']
        logger.info("creating windows")
        regions = targets.window_maker(w=w, g=input.chrom_sizes).intersect(targets)
        p = config['pybedtools']['make_bed_windows']['partitions']
        try:
            assert len(regions) >= p
        except AssertionError as e:
            logger.warning("Number of regions smaller than number of partitions: '{} < {}': lower the number of partitions (config['pybedtools']['make_bed_windows']['partitions']) ".format(len(regions), p))
            raise
        ll = [len(x) for x in regions]
        logger.info("sorting regions")
        ix = sorted(range(len(regions)), key=lambda k: len(regions[k]), reverse=True)
        out = [[regions[i]] for i in ix[0:p]]
        logger.info("grouping regions")
        for j in ix[p:len(ix)]:
            imin = sorted(range(len(out)), key=lambda k: sum(len(r) for r in out[k]))[0]
            out[imin].append(regions[j])
        logger.info("printing regions")
        for j in range(p):
            logger.info("  region {}".format(j))
            b = pybedtools.BedTool("\n".join(str(r) for r in out[j]), from_string=True)
            bedout = pybedtools.BedTool("\n".join(str(r) for r in out[j]), from_string=True)
            bedout.saveas(output.targets[j])
