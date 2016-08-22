#!/bin/bash
export VCF=~/Data/Collaborations/HIVEC/vcf/exome_calls_anno.vcf.gz
export NORMVCF=~/Data/Collaborations/HIVEC/vcf/normalized_exome_calls_anno.vcf.gz
export REF=~/Data/References/Genomes/Homo_sapiens/hg19/Ensembl/Homo_sapiens.GRCh37.75.dna.toplevel.fa
export SNPEFFJAR=/usr/local/bin/snpeff
export db=hiv_exome_seq

zless $VCF \
   | sed 's/ID=AD,Number=./ID=AD,Number=R/' \
   | vt decompose -s - \
   | vt normalize -r $REF - \
   | $SNPEFFJAR GRCh37.75 -formatEff -classic \
   | bgzip -c > $NORMVCF
tabix -p vcf $NORMVCF

gemini load --cores 7 -t snpEff -v $NORMVCF $db

#!/bin/bash
export VCF=~/Data/Collaborations/HIVEC/R_Analysis/exome_calls_anno_snvs_filtered.vcf.bgz
export NORMVCF=~/Data/Collaborations/HIVEC/R_Analysis/normalized_exome_calls_anno_snvs_filtered.vcf.gz
export REF=~/Data/References/Genomes/Homo_sapiens/hg19/Ensembl/Homo_sapiens.GRCh37.75.dna.toplevel.fa
export SNPEFFJAR=/usr/local/bin/snpeff
export db=hiv_exome_seq_snvs_filtered.db

zless $VCF \
   | sed 's/ID=AD,Number=./ID=AD,Number=R/' \
   | vt decompose -s - \
   | vt normalize -r $REF - \
   | $SNPEFFJAR GRCh37.75 -formatEff -classic \
   | bgzip -c > $NORMVCF
tabix -p vcf $NORMVCF

gemini load --cores 7 -t snpEff -v $NORMVCF $db


#!/bin/bash
export VCF=~/Data/Collaborations/HIVEC/R_Analysis/exome_calls_anno_indels_filtered.vcf.bgz
export NORMVCF=~/Data/Collaborations/HIVEC/R_Analysis/normalized_exome_calls_anno_indels_filtered.vcf.gz
export REF=~/Data/References/Genomes/Homo_sapiens/hg19/Ensembl/Homo_sapiens.GRCh37.75.dna.toplevel.fa
export SNPEFFJAR=/usr/local/bin/snpeff
export db=hiv_exome_seq_indels_filtered.db

zless $VCF \
   | sed 's/ID=AD,Number=./ID=AD,Number=R/' \
   | vt decompose -s - \
   | vt normalize -r $REF - \
   | $SNPEFFJAR GRCh37.75 -formatEff -classic \
   | bgzip -c > $NORMVCF
tabix -p vcf $NORMVCF

gemini load --cores 8 -t snpEff -v $NORMVCF $db
