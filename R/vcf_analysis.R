require(VariantAnnotation)
require(TxDb.Hsapiens.UCSC.hg19.knownGene)
require(org.Hs.eg.db)
require(gdata)
library(ensemblVEP)

setwd("~/Data/Collaborations/HIVEC/")
exomes.file <- "./vcf/exome_calls_anno.vcf.gz"
exomes.hdr <- scanVcfHeader(exomes.file)
exomes.hdr.info <- data.frame(info(exomes.hdr))


gr1 <- GRanges(7, IRanges(1, 159138663), strand = "*")

param <- ScanVcfParam(which = gr1, info = c("AC", "AF", "AN", "CCC", "1000g2014oct_all", "ExAC_ALL"), geno = c("GT", "GQ"))
exomes.vcf <- readVcf(exomes.file, "hg19", param = param)
exomes.geno <- readGeno(exomes.file, "GT", param = param)
exomes.info <- readInfo(exomes.file, "AF", param = param)

# post processing of Ensembl Variant Effect Predictor Output (VEP) --------
vep_output_file <- "./vcf/variant_effect_output.txt"
vep <- read.table(vep_output_file, sep = "\t", header = F, stringsAsFactors = F)
colnames(vep) <- c("Uploaded_variation", 
                   "Location", 
                   "Allele", 
                   "Gene", 
                   "Feature", 
                   "Feature_type", 
                   "Consequence",
                   "cDNA_position",
                   "CDS_position",
                   "Protein_position",
                   "Amino_acids",
                   "Codons",
                   "Existing_variation",
                   "Extra")


extrasVector <- c(IMPACT = NA, 
                  DISTANCE = NA,
                  STRAND = NA,
                  VARIANT_CLASS = NA,
                  SYMBOL = NA,
                  SYMBOL_SOURCE = NA,
                  HGNC_ID = NA,
                  BIOTYPE = NA,
                  CANONICAL = NA,
                  TSL = NA,
                  APPRIS = NA,
                  CCDS = NA,
                  ENSP = NA,
                  SWISSPROT = NA,
                  TREMBL = NA,
                  UNIPARC = NA,
                  GENE_PHENO = NA,
                  SIFT = NA,
                  PolyPhen = NA,
                  EXON = NA,
                  INTRON = NA,
                  DOMAINS = NA,
                  GMAF = NA,
                  AFR_MAF = NA,
                  AMR_MAF = NA,
                  EAS_MAF = NA,
                  EUR_MAF = NA,
                  SAS_MAF = NA,
                  AA_MAF = NA,
                  EA_MAF = NA,
                  ExAC_MAF = NA,
                  ExAC_Adj_MAF = NA,
                  ExAC_AFR_MAF = NA,
                  ExAC_AMR_MAF = NA,
                  ExAC_EAS_MAF = NA,
                  ExAC_FIN_MAF = NA,
                  ExAC_NFE_MAF = NA,
                  ExAC_OTH_MAF = NA,
                  ExAC_SAS_MAF = NA,
                  CLIN_SIG = NA,
                  SOMATIC = NA,
                  PHENO = NA,
                  PUBMED = NA,
                  MOTIF_NAME = NA,
                  MOTIF_POS = NA,
                  HIGH_INF_POS = NA,
                  MOTIF_SCORE_CHANGE = NA)

parseExtras <- function(x, y) {
  xv <- extrasVector
  l1 <- length(unlist(strsplit(unlist(strsplit(vep[1,14], ";")), "=")))
  keys <- unlist(strsplit(unlist(strsplit(x, ";")), "="))[seq(1, l1-1, 2)]
  vals <- unlist(strsplit(unlist(strsplit(x, ";")), "="))[seq(2, l1, 2)]
  xv[keys] <- vals
  return(xv)
  }

sapply(seq_along(1:10), function(x){
  v1 <- parseExtras(vep[x, 14])
  length(v1)
  if (x == 1){
    tab1 <- data.frame(matrix(nrow = 1, ncol = length(v1)), stringsAsFactors = F)
    colnames(tab1) <- names(v1)
    tab1[names(v1)] <- v1
  } else {
    tab1[x, names(v1)] <- v1
  }
  return(tab1)
})




