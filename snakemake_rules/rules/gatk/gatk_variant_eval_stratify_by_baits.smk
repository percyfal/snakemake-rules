# -*- snakemake -*-
include: 'gatk.settings'
include: "gatk_variant_eval.rule"

rule gatk_variant_eval_stratify_by_baits:
    """Run GATK VariantEval on a region and stratify results by baits

    The input should be a vcf on which select variants has been run to
    select variants in a given region, typically a gene with exons and
    introns. The stratification regions correspond to baits used
    in a sequence capture.

    """
    wildcard_constraints:
        suffix = "(.vcf|.vcf.gz)"
    params: cmd = config['gatk']['cmd'] + " -T " + VARIANT_EVAL,
            options = " ".join(["-R", config['gatk']['variant_eval']['ref'],
                                config['gatk']['variant_eval']['options'],
                                "--dbsnp {known}".format(known=config['gatk']['known_sites'] if not config['gatk']['known_sites'] == "" else "")]),
            runtime = config['gatk']['variant_eval']['runtime'],
    input: vcf="{prefix}.region_{region}{suffix}", bed="{prefix}.region_{region}.baits.bed"
    output: "{prefix}.region_{region}{suffix}eval_metrics"
    threads: config['gatk']['variant_eval']['threads']
    conda: "env.yaml"
    shell: "{params.cmd} {params.options} --eval {input.vcf} -o {output} --stratIntervals {input.bed} -ST IntervalStratification"

ruleorder: gatk_variant_eval_stratify_by_baits > gatk_variant_eval

