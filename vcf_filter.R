require(VariantTools)
require(VariantAnnotation)
require(GenomicFeatures)

gene.models <- makeTxDbFromGFF("~/Data/References/Annotations/Homo_sapiens/hg19/Ensembl/Homo_sapiens.GRCh37.73.gtf")

prefilters <- FilterRules(list(onlyVariants=function(text){
  !grepl("0/0", text, fixed = TRUE)
}))

vcf.file <- "/Users/u1001407/Data/Collaborations/HIVEC/vcf/exome_calls_anno.vcf.gz"
head(isIndel(vcf))

filtersSNVs <- FilterRules(list(onlySNVs=isSNV))
filtersIndels <- FilterRules(list(onlyIndels=isIndel))

prefilters <- FilterRules(list(onlyVariants=function(text){
!grepl("0/0", text, fixed = TRUE)
}))

filterVcf(vcf.file, 
          genome = "hg19",
          "exome_calls_anno_snvs_filtered.vcf",
          index = TRUE,
          prefilters = prefilters,
          filters = filtersSNVs,
          param=ScanVcfParam(info=NA))

filterVcf(vcf.file, 
          genome = "hg19",
          "exome_calls_anno_indels_filtered.vcf",
          index = TRUE,
          prefilters = prefilters,
          filters = filtersIndels,
          param=ScanVcfParam(info=NA))

indels <- readVcf("exome_calls_anno_indels_filtered.vcf.bgz", genome = "hg19")
snvs <- readVcf("exome_calls_anno_snvs_filtered.vcf.bgz", genome = "hg19")
writeVcf(snvs, "exome_calls_anno_snvs_filtered.vcf", index = FALSE)
writeVcf(snvs, "exome_calls_anno_indels_filtered.vcf", index = FALSE)

testVcf <-readVcf("exome_calls_anno_snvs_filtered.vcf", genome = "hg19")
print(object.size(vcf), unit = "auto")
table(unlist(strsplit(filt(snvs), ";")))

locations <- locateVariants(snvs, gene.models, CodingVariants())
gene_ids <- sub("GeneID:", "",
                locations$GENEID[!is.na(locations$GENEID)])
entrezIDs <- unlist(mget(gene_ids, org.Hs.egENSEMBL2EG, ifnotfound = NA))
entrezIDs <- entrezIDs[!is.na(entrezIDs)]
syms <- unlist(mget(as.character(entrezIDs),
                    org.Hs.egSYMBOL,
                    ifnotfound=NA))

ensemblHost <- "mar2016.archive.ensembl.org"
human <- biomaRt::useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", host = ensemblHost)
ensGenes <- biomaRt::getBM(attributes = c("ensembl_gene_id",
                                          "external_gene_name",
                                          "chromosome_name",
                                          "start_position",
                                          "end_position",
                                          "strand",
                                          "band",
                                          "description",
                                          "percentage_gc_content",
                                          "gene_biotype"),
                           mart = human,
                           filter = "ensembl_gene_id",
                           values = locations$GENEID)
rownames(ensGenes) <- ensGenes$ensembl_gene_id
locations$SYMBOL <- ensGenes[locations$GENEID,]$external_gene_name
locations

