# -*- snakemake -*-
# Simplest greedy algorithm: see https://en.wikipedia.org/wiki/Partition_problem#The_greedy_algorithm
# 1. choose the p largest regions to initial partitions
# 2. for each remaining element, add it to the smallest set
include: "pybedtools.settings.smk"

config_default = {
    'pybedtools' : {
        'make_bed_targets' : {
            'partitions' : 5,
        },
    },
}

update_config(config_default, config)
config = config_default


rule pybedtools_make_bed_targets:
    """Generate bed targets"""
    wildcard_constraints:
        partition = "\d+"
    params: runtime = "01:00:00"
    input: targets = "{prefix}.bed"
    output: targets = ["{{prefix}}.{partition}.bed".format(partition=p+1) for p in range(config['pybedtools']['make_bed_targets']['partitions'])]
    threads: 1
    run:
        regions = pybedtools.BedTool(input.targets)
        p = config['pybedtools']['make_bed_targets']['partitions']
        try:
            assert len(regions) >= p
        except AssertionError as e:
            logger.warning("Number of regions smaller than number of partitions: '{} < {}': lower the number of partitions (config['pybedtools']['make_bed_targets']['partitions']) ".format(len(regions), p))
            raise
        ll = [len(x) for x in regions]
        ix = sorted(range(len(regions)), key=lambda k: len(regions[k]), reverse=True)
        out = [[regions[i]] for i in ix[0:p]]
        for j in ix[p:len(ix)]:
            imin = sorted(range(len(out)), key=lambda k: sum(len(r) for r in out[k]))[0]
            out[imin].append(regions[j])
        for j in range(p):
            bedout = pybedtools.BedTool("\n".join(str(r) for r in out[j]), from_string=True)
            bedout.saveas(output.targets[j])
