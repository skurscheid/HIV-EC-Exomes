require(VariantAnnotation)
require(TxDb.Hsapiens.UCSC.hg19.knownGene)
require(org.Hs.eg.db)
require(gdata)
library(ensemblVEP)

setwd("~/Data/Collaborations/HIVEC/")
exomes.file <- "./vcf/exome_calls_anno.vcf.gz"
exomes.hdr <- scanVcfHeader(exomes.file)

genesym <- c("TRPV1", "TRPV2", "TRPV3", "CXCR4", "CXCR5", "IDH1", "TP53")
genesym <- "CXCR4"
geneid <- select(org.Hs.eg.db, keys=genesym, keytype="SYMBOL",
columns="ENTREZID")

txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
txdb <- renameSeqlevels(txdb, gsub("chr", "", seqlevels(txdb)))

txbygene <- transcriptsBy(txdb, "gene")

gnrng <- unlist(range(txbygene[geneid$ENTREZID]), use.names=FALSE)
names(gnrng) <- geneid$SYMBOL

param <- ScanVcfParam(which = gnrng, info = "DP", geno = c("GT", "GQ"))
vcf <- readVcf(exomes.file, "hg19", param)

dest <- tempfile()
writeVcf(vcf, dest)
VEPparams <- VEPParam()
scriptPath(VEPparams) <- "/usr/local/bin/variant_effect_predictor"
gr <- ensemblVEP(file = dest, param = VEPparams)


# reading in Jim Auer’s analysis results ----------------------------------
r5_resistant_gene_burden <- read.xls("tables/R5_resistant_gene_burden_analysis_fisher.xlsx")

