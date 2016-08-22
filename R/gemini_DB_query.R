require(RSQLite)
require(gdata)

setwd("~/Data/Collaborations/HIVEC/R_Analysis/")

DBs <- c("../geminiDB/hiv_exome_seq", "hiv_exome_seq_snvs_filtered.db", "hiv_exome_seq_indels_filtered.db")
con <- dbConnect(SQLite(), "../geminiDB/hiv_exome_seq")
conFilteredSNVs <- dbConnect(SQLite(), "hiv_exome_seq_snvs_filtered.db")
conFilteredInDels <- dbConnect(SQLite(), "hiv_exome_seq_indels_filtered.db")

con <- conFilteredSNVs

samples <- dbGetQuery(con, "SELECT * from samples")
ped.xls <- read.xls("../geminiDB/PED_table.xlsx", sep = "\t")x
colnames(ped.xls)
colnames(smp)
head(smp)

r5Resistant <- read.table("../ethnicity/R5-resistant.txt", header = F, as.is = T, sep ="\t")
merged <- merge(smp[, c(1,2,3,4,5)], ped.xls[,c(2,5,6,7,8)], by.x = "name", by.y = "name", all.x = T)
merged[is.na(merged$phenotype), ]$phenotype <- 1
merged[is.na(merged$sex), ]$sex <- 3
merged$sex <- as.integer(merged$sex)
merged$hiv_status <- as.character(merged$hiv_status)
merged$ethnicity <- as.character(merged$ethnicity)
merged[is.na(merged$hiv_status), ]$hiv_status <- "none"
merged[is.na(merged$ethnicity), ]$ethnicity <- "none"
merged[merged$ethnicity == -9, ]$ethnicity <- "none"
merged[merged$name %in% r5Resistant$V1, ]$phenotype <- 2
merged$phenotype <- as.integer(merged$phenotype)

write.table(merged, file = "../geminiDB/PED_table.txt", quote = F, sep = "\t", row.names = F, col.names = T)

sapply(DBs, function(x){
  con <- dbConnect(SQLite(), x)
  dbBegin(con)
  dbSendPreparedQuery(con, "update samples set family_id = ?, name = ?, paternal_id = ?, maternal_id = ?, sex = ?, phenotype = ?, hiv_status = ?, ethnicity = ? where sample_id = ?", 
                      bind.data = data.frame(merged$family_id, merged$name, merged$paternal_id, merged$maternal_id, merged$sex, merged$phenotype, merged$hiv_status, merged$ethnicity, merged$sample_id))
  dbCommit(con)
})

samples <- dbGetQuery(con, "SELECT * from samples")

res <- dbGetQuery(con, "select filter, num_hom_ref, num_het, num_hom_alt from variants")
filteredResults <- lapply(DBs[2:3], function(x){
  con <- dbConnect(SQLite(), x)
  res <- dbGetQuery(con,
                    "select chrom,
                    start,
                    end, 
                    variant_id,
                    vcf_id,
                    qual ,
                    type ,
                    sub_type ,
                    call_rate ,
                    in_dbsnp ,
                    rs_ids,
                    in_omim ,
                    clinvar_sig ,
                    clinvar_disease_name ,
                    clinvar_causal_allele ,
                    clinvar_gene_phenotype ,
                    geno2mp_hpo_ct ,
                    pfam_domain ,
                    cyto_band ,
                    rmsk text ,
                    in_cpg_island ,
                    in_segdup ,
                    is_conserved ,
                    num_hom_ref ,
                    num_het ,
                    num_hom_alt ,
                    num_unknown ,
                    aaf ,
                    hwe ,
                    recomb_rate ,
                    gene ,
                    transcript ,
                    is_exonic ,
                    is_coding ,
                    is_splicing ,
                    is_lof ,
                    exon ,
                    codon_change ,
                    aa_change ,
                    aa_length ,
                    biotype,
                    impact,
                    impact_so ,
                    impact_severity ,
                    in_hm2 ,
                    in_hm3 ,
                    is_somatic ,
                    somatic_score ,
                    in_esp ,
                    CAST (aaf_esp_ea AS float) AS aaf_esp_ea,
                    CAST (aaf_esp_aa AS float) AS aaf_esp_aa,
                    CAST (aaf_esp_all AS float) AS aaf_esp_all,
                    exome_chip,
                    in_1kg,
                    CAST (aaf_1kg_amr AS float) AS aaf_1kg_amr,
                    CAST (aaf_1kg_eas AS float) AS aaf_1kg_eas,
                    CAST (aaf_1kg_sas AS float) AS aaf_1kg_sas,
                    CAST (aaf_1kg_afr AS float) AS aaf_1kg_afr,
                    CAST (aaf_1kg_eur AS float) AS aaf_1kg_eur,
                    CAST (aaf_1kg_all AS float) AS aaf_1kg_all,
                    in_exac,
                    CAST (aaf_exac_all AS float) AS aaf_exac_all,
                    CAST (aaf_adj_exac_all AS float) AS aaf_adj_exac_all,
                    CAST (aaf_adj_exac_afr AS float) AS aaf_adj_exac_afr,
                    CAST (aaf_adj_exac_amr AS float) AS aaf_adj_exac_amr,
                    CAST (aaf_adj_exac_eas AS float) AS aaf_adj_exac_eas,
                    CAST (aaf_adj_exac_fin AS float) AS aaf_adj_exac_fin,
                    CAST (aaf_adj_exac_nfe AS float) AS aaf_adj_exac_nfe,
                    CAST (aaf_adj_exac_oth AS float) AS aaf_adj_exac_oth,
                    CAST (aaf_adj_exac_sas AS float) AS aaf_adj_exac_sas,
                    CAST (exac_num_het AS float) AS exac_num_het,
                    CAST (exac_num_hom_alt AS float) AS exac_num_hom_alt,
                    CAST (exac_num_chroms AS float) AS exac_num_chroms, 
                    CAST (max_aaf_all AS float) AS max_aaf_all 
                    from variants
                    where call_rate > 0.99 
                    and impact != 'synonymous_variant'
                    and impact != 'intergenic_variant'
                    and impact != 'intron_variant'
                    ")
})
#where impact_severity != 'LOW' and in_dbsnp = 0 

names(filteredResults) <- DBs[2:3]
filteredVariants <- do.call("rbind", filteredResults)


snv_carriers <- read.table("snv_carriers.txt", header = T, as.is = T, sep = "\t")
