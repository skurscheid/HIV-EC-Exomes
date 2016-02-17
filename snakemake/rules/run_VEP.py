__author__ = "Sebastian Kurscheid (sebastian.kurscheid@anu.edu.au)"
__license__ = "MIT"
__date__ = "2016-02-17"

# vim: syntax=python tabstop=4 expandtab
# coding: utf-8


"""
Rules for running Ensembl variant_effect_predictor.

For usage, include this in your workflow.
"""

configfile: "/Users/u1001407/Development/Collaborations/HIV-EC-Exomes/snakemake/config/config.json"

rule all:
    input:
        expand("variant_effect_predictor/{sample}.vep_output.txt", sample = config["samples"])

rule run_vep:
    input:
        "vcf/exome_calls_anno.{sample}.vcf.gz"
    output:
        "variant_effect_predictor/{sample}.vep_output.txt"
    shell:
        """
            /Users/u1001407/perl5/perlbrew/perls/5.14.4/bin/perl \
            /Users/u1001407/Bioinformatics/ensembl-tools-release-83/scripts/variant_effect_predictor/variant_effect_predictor.pl --input_file  {input} \
                                                                                                                                 --cache \
                                                                                                                                 --offline \
                                                                                                                                 --force_overwrite \
                                                                                                                                 --check_existing \
                                                                                                                                 --check_alleles \
                                                                                                                                 --everything \
                                                                                                                                 --output_file {output}

        """
