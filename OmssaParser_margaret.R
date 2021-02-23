args=(commandArgs(TRUE))

# params
xmlFile <- args[1]
cat(system("date"))
cat(xmlFile)
cat("\n")
# xmlFile <- "/home/bioinformatica/datos/03_Analysis/agarin/13_PME11_May16/05_4SE_2018/Others/L17_PME11_A1_SMIM_Enrich"
require(XML)
data <- xmlParse(paste0(xmlFile, ".xml"))
xml_data_raw <- xmlToList(data)
save(xml_data_raw, file = paste(xmlFile, "_raw.rda", sep = ""))
xml_data_unlist <- unlist(xml_data_raw)
# write.table(xml_data_df, file = "/home/bioinformatica/datos/03_Analysis/agarin/13_PME11_May16/05_4SE_2018/xml_data_df.txt", quote = FALSE, row.names = FALSE, col.names = TRUE, sep="\t")
indexes <- which(names(xml_data_unlist) == names(xml_data_unlist)[1])

xml_data_df <- data.frame(number = numeric(0), evalue = numeric(0), pvalue = numeric(0), charge = numeric(0), peptide = character(0), mass = numeric(0), modType = numeric(0), modSite = numeric(0), theomass = numeric(0), ids = character(0), settingid = numeric(0), MSPepHit_start = numeric(0), MSPepHit_stop = numeric(0), MSPepHit_accession = character(0), MSPepHit_defline = character(0), MSPepHit_protlength = numeric(0), MSPepHit_oid = numeric(0))

for (i in 1:length(indexes)) {
	index <- indexes[i]
	if (i == length(indexes)) {
			index_fin <- which(names(xml_data_unlist) == "MSResponse_scale")
		} else {
			index_fin <- indexes[i+1] - 1
		}
	tmp <- xml_data_unlist[index:index_fin]
	tmp_df <- data.frame(name = names(tmp), value = tmp)
	repetir <- length(tmp_df[tmp_df$name == "MSResponse_hitsets.MSHitSet.MSHitSet_hits.MSHits.MSHits_evalue","value"])
	indexes2 <- which(tmp_df$name == "MSResponse_hitsets.MSHitSet.MSHitSet_hits.MSHits.MSHits_evalue")
	numbervalue <- as.numeric(paste(tmp_df[tmp_df$name == "MSResponse_hitsets.MSHitSet.MSHitSet_number","value"]))
	idsvalue <- paste(tmp_df[tmp_df$name == "MSResponse_hitsets.MSHitSet.MSHitSet_ids.MSHitSet_ids_E","value"])
	settingvalue <- as.numeric(paste(tmp_df[tmp_df$name == "MSResponse_hitsets.MSHitSet.MSHitSet_settingid","value"]))
	if(!(identical(indexes2, integer(0)))) {
		for (j in 1:repetir) {
			index2 <- indexes2[j]
			if (j == length(indexes2)) {
				index2_fin <- which(tmp_df$name == "MSResponse_hitsets.MSHitSet.MSHitSet_hits.MSHits.MSHits_theomass")
				index2_fin <- index2_fin[which(index2_fin > index2)]
			} else {
				index2_fin <- indexes2[j+1] - 1
			}
			tmp2_df <- tmp_df[index2:index2_fin,]
			repetir2 <- length(tmp2_df[tmp2_df$name == "MSResponse_hitsets.MSHitSet.MSHitSet_hits.MSHits.MSHits_pephits.MSPepHit.MSPepHit_start","value"])
			if(repetir2 != 0) {
				out_df_tmp <- data.frame(number = rep(numbervalue, repetir2), evalue = rep(as.numeric(paste(tmp2_df[tmp2_df$name == "MSResponse_hitsets.MSHitSet.MSHitSet_hits.MSHits.MSHits_evalue","value"])), repetir2), pvalue = rep(as.numeric(paste(tmp2_df[tmp2_df$name == "MSResponse_hitsets.MSHitSet.MSHitSet_hits.MSHits.MSHits_pvalue","value"])), repetir2), charge = rep(as.numeric(paste(tmp2_df[tmp2_df$name == "MSResponse_hitsets.MSHitSet.MSHitSet_hits.MSHits.MSHits_charge","value"])), repetir2), peptide = rep(paste(tmp2_df[tmp2_df$name == "MSResponse_hitsets.MSHitSet.MSHitSet_hits.MSHits.MSHits_pepstring","value"]), repetir2), mass = rep(as.numeric(paste(tmp2_df[tmp2_df$name == "MSResponse_hitsets.MSHitSet.MSHitSet_hits.MSHits.MSHits_mass","value"])), repetir2), theomass = rep(as.numeric(paste(tmp2_df[tmp2_df$name == "MSResponse_hitsets.MSHitSet.MSHitSet_hits.MSHits.MSHits_theomass","value"])), repetir2), ids = rep(idsvalue, repetir2), settingid = rep(settingvalue, repetir2), MSPepHit_start = as.numeric(paste(tmp2_df[tmp2_df$name == "MSResponse_hitsets.MSHitSet.MSHitSet_hits.MSHits.MSHits_pephits.MSPepHit.MSPepHit_start","value"])), MSPepHit_stop = as.numeric(paste(tmp2_df[tmp2_df$name == "MSResponse_hitsets.MSHitSet.MSHitSet_hits.MSHits.MSHits_pephits.MSPepHit.MSPepHit_stop","value"])), MSPepHit_accession = paste(tmp2_df[tmp2_df$name == "MSResponse_hitsets.MSHitSet.MSHitSet_hits.MSHits.MSHits_pephits.MSPepHit.MSPepHit_accession","value"]), MSPepHit_defline = paste(tmp2_df[tmp2_df$name == "MSResponse_hitsets.MSHitSet.MSHitSet_hits.MSHits.MSHits_pephits.MSPepHit.MSPepHit_defline","value"]), MSPepHit_protlength = as.numeric(paste(tmp2_df[tmp2_df$name == "MSResponse_hitsets.MSHitSet.MSHitSet_hits.MSHits.MSHits_pephits.MSPepHit.MSPepHit_protlength","value"])), MSPepHit_oid = as.numeric(paste(tmp2_df[tmp2_df$name == "MSResponse_hitsets.MSHitSet.MSHitSet_hits.MSHits.MSHits_pephits.MSPepHit.MSPepHit_oid","value"])))
				repetir3 <- length(tmp2_df[tmp2_df$name == "MSResponse_hitsets.MSHitSet.MSHitSet_hits.MSHits.MSHits_mods.MSModHit.MSModHit_modtype.MSMod","value"])
				if (repetir3 == 0) {
					out_df <- cbind(out_df_tmp, data.frame(modType = rep(-1, repetir2), modSite = rep(-1, repetir2)))
					#cat("option 1. repetir3 0. i=", i, " - j=", j, "\n")
				} else if (repetir3 == 1) {
					out_df <- cbind(out_df_tmp, data.frame(modType = rep(as.numeric(paste(tmp2_df[tmp2_df$name == "MSResponse_hitsets.MSHitSet.MSHitSet_hits.MSHits.MSHits_mods.MSModHit.MSModHit_modtype.MSMod","value"])), repetir2), modSite = rep(as.numeric(paste(tmp2_df[tmp2_df$name == "MSResponse_hitsets.MSHitSet.MSHitSet_hits.MSHits.MSHits_mods.MSModHit.MSModHit_site","value"])), repetir2)))
					#cat("option 2. repetir3 1. i=", i, " - j=", j, "\n")
				} else {
					out_df_tmp2 <- data.frame(modType = as.numeric(paste(tmp2_df[tmp2_df$name == "MSResponse_hitsets.MSHitSet.MSHitSet_hits.MSHits.MSHits_mods.MSModHit.MSModHit_modtype.MSMod","value"])), modSite = as.numeric(paste(tmp2_df[tmp2_df$name == "MSResponse_hitsets.MSHitSet.MSHitSet_hits.MSHits.MSHits_mods.MSModHit.MSModHit_site","value"])))
					out_df <- merge(out_df_tmp, out_df_tmp2)
					#cat("option 3. repetir3>1. i=", i, " - j=", j, "\n")
				}
				xml_data_df <- rbind(xml_data_df, out_df)
			} else {
				cat("careful. there is no start position")
			}
			rm(tmp2_df)
		}
	}
	rm(tmp_df)
}

xml_data_df$modName <- "NA"
xml_data_df[as.numeric(paste(xml_data_df$modType)) == 1,"modName"] <- "Oxidation of M"
xml_data_df[as.numeric(paste(xml_data_df$modType)) == 6,"modName"] <- "Phosphorylation of S"
xml_data_df[as.numeric(paste(xml_data_df$modType)) == 7,"modName"] <- "Phosphorylation of T"
xml_data_df[as.numeric(paste(xml_data_df$modType)) == 8,"modName"] <- "Phosphorylation of Y"
xml_data_df[as.numeric(paste(xml_data_df$modType)) == 10,"modName"] <- "Acetylation of protein N-term"
xml_data_df[as.numeric(paste(xml_data_df$modType)) == 3,"modName"] <- "carbamidomethyl C"

save(xml_data_df, file = paste(xmlFile, ".rda", sep = ""))
write.table(xml_data_df, file = paste(xmlFile, ".txt", sep = ""), quote = FALSE, row.names = FALSE, col.names = TRUE, sep="\t")
cat(system("date"))
