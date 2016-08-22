#!/bin/bash
export VCF=~/Data/Collaborations/HIVEC/vcf/exome_calls_anno.vcf.gz
export NORMVCF=~/Data/Collaborations/HIVEC/vcf/normalized_exome_calls_anno.vcf.gz
export REF=~/Data/References/Genomes/Homo_sapiens/GRCh37_hg19_ensembl75/toplevel/Homo_sapiens.GRCh37.75.dna.toplevel.fa
export SNPEFFJAR=~/Bioinformatics/snpEff/snpEff.jar
export SNPEFFDATA=~/Data/Collaborations/snpEff/data
export VT=~/Bioinformatics/vt/vt
export db=hiv_exome_seq

zless $VCF \
   | sed 's/ID=AD,Number=./ID=AD,Number=R/' \
   | $VT decompose -s - \
   | $VT normalize -r $REF - \
   | java -Xms32G -Xmx40G -jar $SNPEFFJAR GRCh37.75 -formatEff -classic -csvStats snpEff_stats.csv -dataDir $SNPEFFDATA\
   | bgzip -c > $NORMVCF
tabix -p vcf $NORMVCF

gemini load --cores 7 -t snpEff -v $NORMVCF $db
