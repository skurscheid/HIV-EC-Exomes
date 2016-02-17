__author__ = "Sebastian Kurscheid (sebastian.kurscheid@anu.edu.au)"
__license__ = "MIT"
__date__ = "2016-02-17"

# vim: syntax=python tabstop=4 expandtab
# coding: utf-8


"""
Rules for making tabix files from binary VCF files

For usage, include this in your workflow.
"""

configfile: "/Users/u1001407/Development/Collaborations/HIV-EC-Exomes/snakemake/config/config.json"

rule all:
    input:
        expand("vcf/exome_calls_anno.{sample}.vcf.gz.tbi", sample = config["samples"])

rule make_tabix:
    input:
        "vcf/exome_calls_anno.{sample}.vcf.gz"
    output:
        "vcf/exome_calls_anno.{sample}.vcf.gz.tbi"
    shell:
        """
            touch {output}; \
            tabix -f --preset vcf {input}
        """
