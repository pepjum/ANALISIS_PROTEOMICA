#!/usr/bin/Rscript

args=(commandArgs(TRUE))
projectDir <- args[1]		# /home/nostromo/data/03_Analysis/agarin/28_PFortes_Shotgun_lncRNA_Feb18/
fastafileroot <- args[2]	# /home/nostromo/data/03_Analysis/agarin/28_PFortes_Shotgun_lncRNA_Feb18/nextProtDB20170801
nextprotfolder <- args[3]	# /home/nostromo/data/03_Analysis/agarin/23_neXtprot_20170801_Nov2017/

library(Biostrings)
library(doBy)
library(reshape2)
library(ggplot2)
library(stringr)
library(UpSetR)

source("/home/margaret/data/01_Rscripts/A_Funciones/funcionesShotgun.R")
source("/home/margaret/data/01_Rscripts/A_Funciones/funcionesVikv2.R")

cat("loading comet file", "\n")
results_comet<-get(load(paste0(projectDir, "Results/results_comet.rda")))
cat("loading mascot file", "\n")
results_mascot<-get(load(paste0(projectDir, "Results/results_mascot.rda")))
cat("loading omssa file", "\n")
results_omssa<-get(load(paste0(projectDir, "Results/results_omssa.rda")))
cat("loading tandem file", "\n")
results_tandem<-get(load(paste0(projectDir, "Results/results_tandem.rda")))

cat("loading database\n")
db_peptidesXProtAll_filterAA <- get(load(paste0(fastafileroot, "_peptidesXProtAll_filterAA.rda")))

tmp_nextprot <- unique(db_peptidesXProtAll_filterAA[, c("Peptide", "ProteinID")])

results_comet$PeptideSeq <- paste(results_comet$PeptideSeq)
results_comet <- unique(results_comet)
# quito cualquier referencia a Swissprot
results_comet$protID <- NULL
results_comet$pepID <- NULL
results_comet$ProteinAccession <- NULL
results_comet <- unique(results_comet)
results_comet_NX <- merge(results_comet,tmp_nextprot, by.x = "PeptideSeq", by.y = "Peptide", all.x = TRUE)
results_comet_NX <- unique(results_comet_NX)

results_mascot$PeptideSeq <- paste(results_mascot$PeptideSeq)
results_mascot <- unique(results_mascot)
# quito cualquier referencia a Swissprot
results_mascot$protID <- NULL
results_mascot$pepID <- NULL
results_mascot$ProteinAccession <- NULL
results_mascot <- unique(results_mascot)
results_mascot_NX <- merge(results_mascot,tmp_nextprot, by.x = "PeptideSeq", by.y = "Peptide", all.x = TRUE)
results_mascot_NX <- unique(results_mascot_NX)

#projectDir<-"/home/nostromo/data/pepe/EMBRIO_5_04/"
results_omssa$PeptideSeq <- paste(results_omssa$PeptideSeq)
results_omssa <- unique(results_omssa)
# quito cualquier referencia a Swissprot
results_omssa$protID <- NULL
results_omssa$pepID <- NULL
results_omssa$ProteinAccession <- NULL
results_omssa$PeptideSeq<-toupper(results_omssa$PeptideSeq)
results_omssa <- unique(results_omssa)
results_omssa_NX <- merge(results_omssa,tmp_nextprot, by.x = "PeptideSeq", by.y = "Peptide", all.x = TRUE)
results_omssa_NX <- unique(results_omssa_NX)

results_tandem$PeptideSeq <- paste(results_tandem$PeptideSeq)
results_tandem <- unique(results_tandem)
# quito cualquier referencia a Swissprot
results_tandem$protID <- NULL
results_tandem$pepID <- NULL
results_tandem$ProteinAccession <- NULL
results_tandem <- unique(results_tandem)
results_tandem_NX <- merge(results_tandem,tmp_nextprot, by.x = "PeptideSeq", by.y = "Peptide", all.x = TRUE)
results_tandem_NX <- unique(results_tandem_NX)

results_comet_NX <- results_comet_NX[!(is.na(results_comet_NX$ProteinID)),]
results_mascot_NX <- results_mascot_NX[!(is.na(results_mascot_NX$ProteinID)),]
results_omssa_NX <- results_omssa_NX[!(is.na(results_omssa_NX$ProteinID)),]
results_tandem_NX <- results_tandem_NX[!(is.na(results_tandem_NX$ProteinID)),]

# results_comet_NX$PSM <- paste(results_comet_NX$datfile, results_comet_NX$Query, sep = "-")
# results_mascot_NX$PSM <- paste(results_mascot_NX$datfile, results_mascot_NX$Query, sep = "-")
# results_omssa_NX$PSM <- paste(results_omssa_NX$datfile, results_omssa_NX$Query, sep = "-")
# results_tandem_NX$PSM <- paste(results_tandem_NX$datfile, results_tandem_NX$Query, sep = "-")


save(results_comet_NX, file = paste0(projectDir, "Results/results_comet_NX.rda"))
save(results_mascot_NX, file = paste0(projectDir, "Results/results_mascot_NX.rda"))
save(results_omssa_NX, file = paste0(projectDir, "Results/results_omssa_NX.rda"))
save(results_tandem_NX, file = paste0(projectDir, "Results/results_tandem_NX.rda"))

# Venn diagrams
results_comet_NX$NextprotID<-sapply(strsplit(paste(results_comet_NX$ProteinID),"-"),"[",1)
results_mascot_NX$NextprotID<-sapply(strsplit(paste(results_mascot_NX$ProteinID),"-"),"[",1)
results_omssa_NX$NextprotID<-sapply(strsplit(paste(results_omssa_NX$ProteinID),"-"),"[",1)
results_tandem_NX$NextprotID<-sapply(strsplit(paste(results_tandem_NX$ProteinID),"-"),"[",1)


pdf(file = paste0(projectDir, "Plots/results_Venn_NXProteinsPerSearchEngine.pdf"), width = 10, height = 10)
compare4List(paste(unique(results_comet_NX$NextprotID)), paste(unique(results_mascot_NX$NextprotID)), paste(unique(results_omssa_NX$NextprotID)), paste(unique(results_tandem_NX$NextprotID)), "Comet results", "Mascot results", "Omssa results", "Tandem results",  "Proteins (Nextprot) for each search engine")
dev.off()

library(UpSetR)

listInput <- list(mascot = paste(results_mascot_NX$NextprotID) , comet = paste(results_comet_NX$NextprotID), omssa =paste(results_omssa_NX$NextprotID), tandem=paste(results_tandem_NX$NextprotID)  )


pdf(file = paste0(projectDir, "Plots/results_upsetR_NextProtID_PerSearchEngine.pdf"), width = 10, height = 10)
upset(fromList(listInput), order.by = "freq")
dev.off()


# discriminant peptides

db_peptidesXProtAll_filterAA_disc <- get(load(paste0(fastafileroot,"_peptidesXProtAll_filterAA_disc.rda")))
load(paste0(nextprotfolder, "nextProtXChrXENSP.RData"))

results_comet_NX_annot <- merge(results_comet_NX, unique(nextProtXChrXENSP[, c("NextprotID", "Missing", "Chr")]), by.x = "NextprotID", by.y = "NextprotID")
results_comet_NX_annot$discriminant <- 0
results_comet_NX_annot[paste(results_comet_NX_annot$PeptideSeq) %in% paste(db_peptidesXProtAll_filterAA_disc$Peptide), "discriminant"] <- 1
results_comet_NX_annot_disc <- results_comet_NX_annot[results_comet_NX_annot$discriminant == 1, ]

results_mascot_NX_annot <- merge(results_mascot_NX, unique(nextProtXChrXENSP[, c("NextprotID", "Missing", "Chr")]), by.x = "NextprotID", by.y = "NextprotID")
results_mascot_NX_annot$discriminant <- 0
results_mascot_NX_annot[paste(results_mascot_NX_annot$PeptideSeq) %in% paste(db_peptidesXProtAll_filterAA_disc$Peptide), "discriminant"] <- 1
results_mascot_NX_annot_disc <- results_mascot_NX_annot[results_mascot_NX_annot$discriminant == 1, ]

results_omssa_NX_annot <- merge(results_omssa_NX, unique(nextProtXChrXENSP[, c("NextprotID", "Missing", "Chr")]), by.x = "NextprotID", by.y = "NextprotID")
results_omssa_NX_annot$discriminant <- 0
results_omssa_NX_annot[paste(results_omssa_NX_annot$PeptideSeq) %in% paste(db_peptidesXProtAll_filterAA_disc$Peptide), "discriminant"] <- 1
results_omssa_NX_annot_disc <- results_omssa_NX_annot[results_omssa_NX_annot$discriminant == 1, ]

results_tandem_NX_annot <- merge(results_tandem_NX, unique(nextProtXChrXENSP[, c("NextprotID", "Missing", "Chr")]), by.x = "NextprotID", by.y = "NextprotID")
results_tandem_NX_annot$discriminant <- 0
results_tandem_NX_annot[paste(results_tandem_NX_annot$PeptideSeq) %in% paste(db_peptidesXProtAll_filterAA_disc$Peptide), "discriminant"] <- 1
results_tandem_NX_annot_disc <- results_tandem_NX_annot[results_tandem_NX_annot$discriminant == 1, ]

save(results_comet_NX_annot, file = paste0(projectDir, "Results/results_comet_NX_annot.rda"))
save(results_comet_NX_annot_disc, file = paste0(projectDir, "Results/results_comet_NX_annot_disc.rda"))
save(results_mascot_NX_annot, file = paste0(projectDir, "Results/results_mascot_NX_annot.rda"))
save(results_mascot_NX_annot_disc, file = paste0(projectDir, "Results/results_mascot_NX_annot_disc.rda"))
save(results_omssa_NX_annot, file = paste0(projectDir, "Results/results_omssa_NX_annot.rda"))
save(results_omssa_NX_annot_disc, file = paste0(projectDir, "Results/results_omssa_NX_annot_disc.rda"))
save(results_tandem_NX_annot, file = paste0(projectDir, "Results/results_tandem_NX_annot.rda"))
save(results_tandem_NX_annot_disc, file = paste0(projectDir, "Results/results_tandem_NX_annot_disc.rda"))

results_comet_NX_annot_disc_missing <- results_comet_NX_annot_disc[results_comet_NX_annot_disc$Missing == 1, ]
results_mascot_NX_annot_disc_missing <- results_mascot_NX_annot_disc[results_mascot_NX_annot_disc$Missing == 1, ]
results_omssa_NX_annot_disc_missing <- results_omssa_NX_annot_disc[results_omssa_NX_annot_disc$Missing == 1, ]
results_tandem_NX_annot_disc_missing <- results_tandem_NX_annot_disc[results_tandem_NX_annot_disc$Missing == 1, ]

pdf(file = paste0(projectDir, "Plots/results_Venn_NX_MissingProteinsPerSearchEngine_before_uniqueness_checker.pdf"), width = 10, height = 10)
compare4List(paste(unique(results_comet_NX_annot_disc_missing$NextprotID)), paste(unique(results_mascot_NX_annot_disc_missing$NextprotID)), paste(unique(results_omssa_NX_annot_disc_missing$NextprotID)), paste(unique(results_tandem_NX_annot_disc_missing$NextprotID)), "Comet results", "Mascot results", "Omssa results", "Tandem results",  "Missing Proteins (Nextprot) for each search engine")
dev.off()

listInput <- list(mascot = paste(results_mascot_NX_annot_disc_missing$NextprotID) , comet = paste(results_comet_NX_annot_disc_missing$NextprotID), omssa =paste(results_omssa_NX_annot_disc_missing$NextprotID), tandem=paste(results_tandem_NX_annot_disc_missing$NextprotID)  )


pdf(file = paste0(projectDir, "Plots/results_upsetR_NextProtID_MissingProteinsPerSearchEngine_before_uniqueness_checker.pdf"), width = 10, height = 10)
upset(fromList(listInput), order.by = "freq")
dev.off()


results_tandem_NX_annot_disc_missing$PSM<-toupper(results_tandem_NX_annot_disc_missing$PSM)
results_comet_NX_annot_disc_missing$PSM<-toupper(results_comet_NX_annot_disc_missing$PSM)
results_omssa_NX_annot_disc_missing$PSM<-toupper(results_omssa_NX_annot_disc_missing$PSM)
results_mascot_NX_annot_disc_missing$PSM<-toupper(results_mascot_NX_annot_disc_missing$PSM)

results_comet_NX_annot_disc_missing$PSM<-gsub("-","_",results_comet_NX_annot_disc_missing$PSM)
results_comet_NX_annot_disc_missing$PSM<-gsub("#","_",results_comet_NX_annot_disc_missing$PSM)

results_mascot_NX_annot_disc_missing$PSM<-gsub("-","_",results_mascot_NX_annot_disc_missing$PSM)
results_mascot_NX_annot_disc_missing$PSM<-gsub("#","_",results_mascot_NX_annot_disc_missing$PSM)

results_omssa_NX_annot_disc_missing$PSM<-gsub("-","_",results_omssa_NX_annot_disc_missing$PSM)
results_omssa_NX_annot_disc_missing$PSM<-gsub("#","_",results_omssa_NX_annot_disc_missing$PSM)

results_tandem_NX_annot_disc_missing$PSM<-gsub("-","_",results_tandem_NX_annot_disc_missing$PSM)
results_tandem_NX_annot_disc_missing$PSM<-gsub("#","_",results_tandem_NX_annot_disc_missing$PSM)



pdf(file = paste0(projectDir, "Plots/results_Venn_NX_MissingProteinsPerSearchEngine_PSMs_before_uniqueness_checker.pdf"), width = 10, height = 10)
compare4List(paste(unique(results_comet_NX_annot_disc_missing$PSM)), paste(unique(results_mascot_NX_annot_disc_missing$PSM)), paste(unique(results_omssa_NX_annot_disc_missing$PSM)), paste(unique(results_tandem_NX_annot_disc_missing$PSM)), "Comet results", "Mascot results", "Omssa results", "Tandem results",  "PSMs of Missing Proteins for each search engine")
dev.off()

listInput <- list(mascot = paste(results_mascot_NX_annot_disc_missing$PSM) , comet = paste(results_comet_NX_annot_disc_missing$PSM), omssa =paste(results_omssa_NX_annot_disc_missing$PSM), tandem=paste(results_tandem_NX_annot_disc_missing$PSM)  )


pdf(file = paste0(projectDir, "Plots/results_upsetR_NX_MissingProteinsPerSearchEngine_PSMs_before_uniqueness_checker.pdf"), width = 10, height = 10)
upset(fromList(listInput), order.by = "freq")
dev.off()



save(results_comet_NX_annot_disc_missing, file = paste0(projectDir, "Results/results_comet_NX_annot_disc_missing.rda"))
save(results_mascot_NX_annot_disc_missing, file = paste0(projectDir, "Results/results_mascot_NX_annot_disc_missing.rda"))
save(results_omssa_NX_annot_disc_missing, file = paste0(projectDir, "Results/results_omssa_NX_annot_disc_missing.rda"))
save(results_tandem_NX_annot_disc_missing, file = paste0(projectDir, "Results/results_tandem_NX_annot_disc_missing.rda"))

write.table(unique(results_comet_NX_annot_disc_missing$PeptideSeq), file = paste0(projectDir, "Results/results_comet_NX_annot_disc_missing_peptides.txt"), col.names=FALSE, row.names = FALSE, sep="\t", quote=FALSE)

write.table(unique(results_mascot_NX_annot_disc_missing$PeptideSeq), file = paste0(projectDir, "Results/results_mascot_NX_annot_disc_missing_peptides.txt"), col.names=FALSE, row.names = FALSE, sep="\t", quote=FALSE)

write.table(unique(results_omssa_NX_annot_disc_missing$PeptideSeq), file = paste0(projectDir, "Results/results_omssa_NX_annot_disc_missing_peptides.txt"), col.names=FALSE, row.names = FALSE, sep="\t", quote=FALSE)

write.table(unique(results_tandem_NX_annot_disc_missing$PeptideSeq), file = paste0(projectDir, "Results/results_tandem_NX_annot_disc_missing_peptides.txt"), col.names=FALSE, row.names = FALSE, sep="\t", quote=FALSE)
