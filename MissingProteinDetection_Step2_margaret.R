#!/usr/bin/Rscript

args=(commandArgs(TRUE))
projectDir <- args[1]		# /home/nostromo/data/03_Analysis/agarin/28_PFortes_Shotgun_lncRNA_Feb18/
fastafileroot <- args[2]	# /home/nostromo/data/03_Analysis/agarin/28_PFortes_Shotgun_lncRNA_Feb18/nextProtDB20170801
nextprotfolder <- args[3]	# /home/nostromo/data/03_Analysis/agarin/23_neXtprot_20170801_Nov17/

library(Biostrings)
library(doBy)
library(reshape2)
library(ggplot2)
library(stringr)
source("/home/margaret/data/01_Rscripts/A_Funciones/funcionesShotgun.R")
source("/home/margaret/data/01_Rscripts/A_Funciones/funcionesVikv2.R")

load(paste0(projectDir, "Results/results_comet_NX_annot_disc_missing.rda"))
load(paste0(projectDir, "Results/results_mascot_NX_annot_disc_missing.rda"))
load(paste0(projectDir, "Results/results_omssa_NX_annot_disc_missing.rda"))
load(paste0(projectDir, "Results/results_tandem_NX_annot_disc_missing.rda"))


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

results_omssa_NX_annot_disc_missing$PeptideSeq<-toupper(results_omssa_NX_annot_disc_missing$PeptideSeq)



results_comet_NX_annot_disc_missing_peptides_uniqueness <- read.csv2(paste0(projectDir, "Results/results_comet_NX_annot_disc_missing_peptides_uniqueness.csv"), header = TRUE, sep = ",")
if(dim(results_comet_NX_annot_disc_missing_peptides_uniqueness)[1] != 0) {
	results_comet_NX_annot_disc_missing_uniqueness <- results_comet_NX_annot_disc_missing[paste(results_comet_NX_annot_disc_missing$PeptideSeq) %in% paste(results_comet_NX_annot_disc_missing_peptides_uniqueness$peptide), ]
	save(results_comet_NX_annot_disc_missing_uniqueness, file = paste0(projectDir, "Results/results_comet_NX_annot_disc_missing_uniqueness.rda"))
	tmp <- unique(results_comet_NX_annot_disc_missing_uniqueness[, c("sample", "NextprotID", "PeptideSeq")])
	#tmp <- unique(results_comet_NX_annot_disc_missing_uniqueness[, c("sample", "ProteinID", "PeptideSeq")])

	tmp2 <- summaryBy(PeptideSeq~NextprotID + sample, data=tmp, FUN=length, keep.names=TRUE)
	#tmp2 <- summaryBy(PeptideSeq~ProteinID + sample, data=tmp, FUN=length, keep.names=TRUE)

	colnames(tmp2)[3] <- "NofDiscPeptides"
	results_comet_NX_annot_disc_missing_uniqueness <- merge(results_comet_NX_annot_disc_missing_uniqueness, tmp2, by.x = c("NextprotID", "sample"), by.y = c("NextprotID", "sample"), all.x = TRUE)
#	results_comet_NX_annot_disc_missing_uniqueness <- merge(results_comet_NX_annot_disc_missing_uniqueness, tmp2, by.x = c("ProteinID", "sample"), by.y = c("ProteinID", "sample"), all.x = TRUE)
	results_comet_NX_annot_disc_missing_uniqueness[is.na(results_comet_NX_annot_disc_missing_uniqueness$NofDiscPeptides), "NofDiscPeptides"] <- 0
	results_comet_NX_annot_disc_missing_uniqueness_one_hit <- results_comet_NX_annot_disc_missing_uniqueness[results_comet_NX_annot_disc_missing_uniqueness$NofDiscPeptides == 1,]

	results_comet_NX_annot_disc_missing_uniqueness_2UniqPep <- results_comet_NX_annot_disc_missing_uniqueness[results_comet_NX_annot_disc_missing_uniqueness$NofDiscPeptides > 1,]
	save(results_comet_NX_annot_disc_missing_uniqueness_2UniqPep, file = paste0(projectDir, "Results/results_comet_NX_annot_disc_missing_uniqueness_2UniqPep.rda"))
	results_comet_NX_annot_disc_missing_uniqueness_f <- unique(results_comet_NX_annot_disc_missing_uniqueness[, c("PeptideSeq", "NextprotID", "search_engine", "Chr", "sample", "NofDiscPeptides", "PSM")])
#	results_comet_NX_annot_disc_missing_uniqueness_f <- unique(results_comet_NX_annot_disc_missing_uniqueness[, c("PeptideSeq", "ProteinID", "search_engine", "Chr", "sample", "NofDiscPeptides")])

	results_comet_NX_annot_disc_missing_uniqueness_ff <- aggregate(sample ~ PeptideSeq + PSM + NextprotID + search_engine + Chr + NofDiscPeptides, data = results_comet_NX_annot_disc_missing_uniqueness_f, FUN = paste, collapse = ",")
#	results_comet_NX_annot_disc_missing_uniqueness_ff <- aggregate(sample ~ PeptideSeq + ProteinID + search_engine + Chr + NofDiscPeptides, data = results_comet_NX_annot_disc_missing_uniqueness_f, FUN = paste, collapse = ",")

} else {
	cat("Comet has not identified any missing proteins!:")
	results_comet_NX_annot_disc_missing_uniqueness <- results_comet_NX_annot_disc_missing[0,]
	save(results_comet_NX_annot_disc_missing_uniqueness, file = paste0(projectDir, "Results/results_comet_NX_annot_disc_missing_uniqueness.rda"))
	results_comet_NX_annot_disc_missing_uniqueness_2UniqPep <- results_comet_NX_annot_disc_missing_uniqueness
	save(results_comet_NX_annot_disc_missing_uniqueness_2UniqPep, file = paste0(projectDir, "Results/results_comet_NX_annot_disc_missing_uniqueness_2UniqPep.rda"))
	results_comet_NX_annot_disc_missing_uniqueness_2UniqPep_f <- results_comet_NX_annot_disc_missing_uniqueness_2UniqPep
	results_comet_NX_annot_disc_missing_uniqueness_2UniqPep_ff <- results_comet_NX_annot_disc_missing_uniqueness_2UniqPep
}

results_mascot_NX_annot_disc_missing_peptides_uniqueness <- read.csv2(paste0(projectDir, "Results/results_mascot_NX_annot_disc_missing_peptides_uniqueness.csv"), header = TRUE, sep = ",")
if(dim(results_mascot_NX_annot_disc_missing_peptides_uniqueness)[1] != 0) {
	results_mascot_NX_annot_disc_missing_uniqueness <- results_mascot_NX_annot_disc_missing[paste(results_mascot_NX_annot_disc_missing$PeptideSeq) %in% paste(results_mascot_NX_annot_disc_missing_peptides_uniqueness$peptide), ]
	save(results_mascot_NX_annot_disc_missing_uniqueness, file = paste0(projectDir, "Results/results_mascot_NX_annot_disc_missing_uniqueness.rda"))
	tmp <- unique(results_mascot_NX_annot_disc_missing_uniqueness[, c("sample", "NextprotID", "PeptideSeq")])
	#tmp <- unique(results_mascot_NX_annot_disc_missing_uniqueness[, c("sample", "ProteinID", "PeptideSeq")])

	tmp2 <- summaryBy(PeptideSeq~NextprotID + sample, data=tmp, FUN=length, keep.names=TRUE)
	#tmp2 <- summaryBy(PeptideSeq~ProteinID + sample, data=tmp, FUN=length, keep.names=TRUE)

	colnames(tmp2)[3] <- "NofDiscPeptides"
	results_mascot_NX_annot_disc_missing_uniqueness <- merge(results_mascot_NX_annot_disc_missing_uniqueness, tmp2, by.x = c("NextprotID", "sample"), by.y = c("NextprotID", "sample"), all.x = TRUE)
#	results_mascot_NX_annot_disc_missing_uniqueness <- merge(results_mascot_NX_annot_disc_missing_uniqueness, tmp2, by.x = c("ProteinID", "sample"), by.y = c("ProteinID", "sample"), all.x = TRUE)

	results_mascot_NX_annot_disc_missing_uniqueness[is.na(results_mascot_NX_annot_disc_missing_uniqueness$NofDiscPeptides), "NofDiscPeptides"] <- 0
	results_mascot_NX_annot_disc_missing_uniqueness_one_hit <- results_mascot_NX_annot_disc_missing_uniqueness[results_mascot_NX_annot_disc_missing_uniqueness$NofDiscPeptides == 1,]

	results_mascot_NX_annot_disc_missing_uniqueness_2UniqPep <- results_mascot_NX_annot_disc_missing_uniqueness[results_mascot_NX_annot_disc_missing_uniqueness$NofDiscPeptides > 1,]
	save(results_mascot_NX_annot_disc_missing_uniqueness_2UniqPep, file = paste0(projectDir, "Results/results_mascot_NX_annot_disc_missing_uniqueness_2UniqPep.rda"))
#	results_mascot_NX_annot_disc_missing_uniqueness_f <- unique(results_mascot_NX_annot_disc_missing_uniqueness[, c("PeptideSeq", "ProteinID", "search_engine", "Chr", "sample", "NofDiscPeptides")])
	results_mascot_NX_annot_disc_missing_uniqueness_f <- unique(results_mascot_NX_annot_disc_missing_uniqueness[, c("PeptideSeq", "NextprotID", "search_engine", "Chr", "sample", "NofDiscPeptides", "PSM")])

#	results_mascot_NX_annot_disc_missing_uniqueness_ff <- aggregate(sample ~ PeptideSeq + ProteinID + search_engine + Chr + NofDiscPeptides, data = results_mascot_NX_annot_disc_missing_uniqueness_f, FUN = paste, collapse = ",")
	results_mascot_NX_annot_disc_missing_uniqueness_ff <- aggregate(sample ~ PeptideSeq + PSM + NextprotID + search_engine + Chr + NofDiscPeptides, data = results_mascot_NX_annot_disc_missing_uniqueness_f, FUN = paste, collapse = ",")

} else  {
	cat("Mascot has not identified any missing proteins!:")
	results_mascot_NX_annot_disc_missing_uniqueness <- results_mascot_NX_annot_disc_missing[0,]
	save(results_mascot_NX_annot_disc_missing_uniqueness, file = paste0(projectDir, "Results/results_mascot_NX_annot_disc_missing_uniqueness.rda"))
	results_mascot_NX_annot_disc_missing_uniqueness_2UniqPep <- results_mascot_NX_annot_disc_missing_uniqueness
	save(results_mascot_NX_annot_disc_missing_uniqueness_2UniqPep, file = paste0(projectDir, "Results/results_mascot_NX_annot_disc_missing_uniqueness_2UniqPep.rda"))
	results_mascot_NX_annot_disc_missing_uniqueness_f <- results_mascot_NX_annot_disc_missing_uniqueness_2UniqPep
	results_mascot_NX_annot_disc_missing_uniqueness_ff <- results_mascot_NX_annot_disc_missing_uniqueness_2UniqPep
}

results_omssa_NX_annot_disc_missing_peptides_uniqueness <- read.csv2(paste0(projectDir, "Results/results_omssa_NX_annot_disc_missing_peptides_uniqueness.csv"), header = TRUE, sep = ",")
if(dim(results_omssa_NX_annot_disc_missing_peptides_uniqueness)[1] != 0) {
	results_omssa_NX_annot_disc_missing_uniqueness <- results_omssa_NX_annot_disc_missing[paste(results_omssa_NX_annot_disc_missing$PeptideSeq) %in% paste(results_omssa_NX_annot_disc_missing_peptides_uniqueness$peptide), ]
	save(results_omssa_NX_annot_disc_missing_uniqueness, file = paste0(projectDir, "Results/results_omssa_NX_annot_disc_missing_uniqueness.rda"))
	tmp <- unique(results_omssa_NX_annot_disc_missing_uniqueness[, c("sample", "NextprotID", "PeptideSeq")])
	#tmp <- unique(results_omssa_NX_annot_disc_missing_uniqueness[, c("sample", "ProteinID", "PeptideSeq")])

	tmp2 <- summaryBy(PeptideSeq~NextprotID + sample, data=tmp, FUN=length, keep.names=TRUE)
	#tmp2 <- summaryBy(PeptideSeq~ProteinID + sample, data=tmp, FUN=length, keep.names=TRUE)

	colnames(tmp2)[3] <- "NofDiscPeptides"
	results_omssa_NX_annot_disc_missing_uniqueness <- merge(results_omssa_NX_annot_disc_missing_uniqueness, tmp2, by.x = c("NextprotID", "sample"), by.y = c("NextprotID", "sample"), all.x = TRUE)
	#results_omssa_NX_annot_disc_missing_uniqueness <- merge(results_omssa_NX_annot_disc_missing_uniqueness, tmp2, by.x = c("ProteinID", "sample"), by.y = c("ProteinID", "sample"), all.x = TRUE)

	results_omssa_NX_annot_disc_missing_uniqueness[is.na(results_omssa_NX_annot_disc_missing_uniqueness$NofDiscPeptides), "NofDiscPeptides"] <- 0
	results_omssa_NX_annot_disc_missing_uniqueness_2UniqPep <- results_omssa_NX_annot_disc_missing_uniqueness[results_omssa_NX_annot_disc_missing_uniqueness$NofDiscPeptides > 1,]
	results_omssa_NX_annot_disc_missing_uniqueness_one_hit <- results_omssa_NX_annot_disc_missing_uniqueness[results_omssa_NX_annot_disc_missing_uniqueness$NofDiscPeptides == 1,]

	save(results_omssa_NX_annot_disc_missing_uniqueness_2UniqPep, file = paste0(projectDir, "Results/results_omssa_NX_annot_disc_missing_uniqueness_2UniqPep.rda"))
	results_omssa_NX_annot_disc_missing_uniqueness_f <- unique(results_omssa_NX_annot_disc_missing_uniqueness[, c("PeptideSeq", "NextprotID", "search_engine", "Chr", "sample", "NofDiscPeptides", "PSM")])
#	results_omssa_NX_annot_disc_missing_uniqueness_f <- unique(results_omssa_NX_annot_disc_missing_uniqueness[, c("PeptideSeq", "ProteinID", "search_engine", "Chr", "sample", "NofDiscPeptides")])

	results_omssa_NX_annot_disc_missing_uniqueness_ff <- aggregate(sample ~ PeptideSeq + PSM + NextprotID + search_engine + Chr + NofDiscPeptides, data = results_omssa_NX_annot_disc_missing_uniqueness_f, FUN = paste, collapse = ",")
#	results_omssa_NX_annot_disc_missing_uniqueness_ff <- aggregate(sample ~ PeptideSeq + ProteinID + search_engine + Chr + NofDiscPeptides, data = results_omssa_NX_annot_disc_missing_uniqueness_f, FUN = paste, collapse = ",")

} else  {
	cat("Omssa has not identified any missing proteins!:")
	results_omssa_NX_annot_disc_missing_uniqueness <- results_omssa_NX_annot_disc_missing[0,]
	save(results_omssa_NX_annot_disc_missing_uniqueness, file = paste0(projectDir, "Results/results_omssa_NX_annot_disc_missing_uniqueness.rda"))
	results_omssa_NX_annot_disc_missing_uniqueness_2UniqPep <- results_omssa_NX_annot_disc_missing_uniqueness

	save(results_omssa_NX_annot_disc_missing_uniqueness_2UniqPep, file = paste0(projectDir, "Results/results_omssa_NX_annot_disc_missing_uniqueness_2UniqPep.rda"))
	results_omssa_NX_annot_disc_missing_uniqueness_f <- results_omssa_NX_annot_disc_missing_uniqueness_2UniqPep
	results_omssa_NX_annot_disc_missing_uniqueness_ff <- results_omssa_NX_annot_disc_missing_uniqueness_2UniqPep
}

results_tandem_NX_annot_disc_missing_peptides_uniqueness <- read.csv2(paste0(projectDir, "Results/results_tandem_NX_annot_disc_missing_peptides_uniqueness.csv"), header = TRUE, sep = ",")
if(dim(results_tandem_NX_annot_disc_missing_peptides_uniqueness)[1] != 0) {
	results_tandem_NX_annot_disc_missing_uniqueness <- results_tandem_NX_annot_disc_missing[paste(results_tandem_NX_annot_disc_missing$PeptideSeq) %in% paste(results_tandem_NX_annot_disc_missing_peptides_uniqueness$peptide), ]
	save(results_tandem_NX_annot_disc_missing_uniqueness, file = paste0(projectDir, "Results/results_tandem_NX_annot_disc_missing_uniqueness.rda"))
	tmp <- unique(results_tandem_NX_annot_disc_missing_uniqueness[, c("sample", "NextprotID", "PeptideSeq")])
	tmp2 <- summaryBy(PeptideSeq~NextprotID + sample, data=tmp, FUN=length, keep.names=TRUE)
#	tmp <- unique(results_omssa_NX_annot_disc_missing_uniqueness[, c("sample", "ProteinID", "PeptideSeq")])
#	tmp2 <- summaryBy(PeptideSeq~ProteinID + sample, data=tmp, FUN=length, keep.names=TRUE)

	colnames(tmp2)[3] <- "NofDiscPeptides"
	results_tandem_NX_annot_disc_missing_uniqueness <- merge(results_tandem_NX_annot_disc_missing_uniqueness, tmp2, by.x = c("NextprotID", "sample"), by.y = c("NextprotID", "sample"), all.x = TRUE)
	#results_tandem_NX_annot_disc_missing_uniqueness <- merge(results_tandem_NX_annot_disc_missing_uniqueness, tmp2, by.x = c("ProteinID", "sample"), by.y = c("ProteinID", "sample"), all.x = TRUE)

	results_tandem_NX_annot_disc_missing_uniqueness[is.na(results_tandem_NX_annot_disc_missing_uniqueness$NofDiscPeptides), "NofDiscPeptides"] <- 0
	results_tandem_NX_annot_disc_missing_uniqueness_2UniqPep <- results_tandem_NX_annot_disc_missing_uniqueness[results_tandem_NX_annot_disc_missing_uniqueness$NofDiscPeptides > 1,]
	results_tandem_NX_annot_disc_missing_uniqueness_one_hit <- results_tandem_NX_annot_disc_missing_uniqueness[results_tandem_NX_annot_disc_missing_uniqueness$NofDiscPeptides == 1,]

	save(results_tandem_NX_annot_disc_missing_uniqueness_2UniqPep, file = paste0(projectDir, "Results/results_tandem_NX_annot_disc_missing_uniqueness_2UniqPep.rda"))
	results_tandem_NX_annot_disc_missing_uniqueness_f <- unique(results_tandem_NX_annot_disc_missing_uniqueness[, c("PeptideSeq", "NextprotID", "search_engine", "Chr", "sample", "NofDiscPeptides","PSM")])
#	results_tandem_NX_annot_disc_missing_uniqueness_f <- unique(results_tandem_NX_annot_disc_missing_uniqueness[, c("PeptideSeq", "ProteinID", "search_engine", "Chr", "sample", "NofDiscPeptides")])
	results_tandem_NX_annot_disc_missing_uniqueness_ff <- aggregate(sample ~ PeptideSeq + PSM + NextprotID + search_engine + Chr + NofDiscPeptides, data = results_tandem_NX_annot_disc_missing_uniqueness_f, FUN = paste, collapse = ",")
#	results_tandem_NX_annot_disc_missing_uniqueness_ff <- aggregate(sample ~ PeptideSeq + ProteinID + search_engine + Chr + NofDiscPeptides, data = results_tandem_NX_annot_disc_missing_uniqueness_f, FUN = paste, collapse = ",")

} else  {
	cat("Tandem has not identified any missing proteins!:")
	results_tandem_NX_annot_disc_missing_uniqueness <- results_tandem_NX_annot_disc_missing[0,]
	save(results_tandem_NX_annot_disc_missing_uniqueness, file = paste0(projectDir, "Results/results_tandem_NX_annot_disc_missing_uniqueness.rda"))
	results_tandem_NX_annot_disc_missing_uniqueness_2UniqPep <- results_tandem_NX_annot_disc_missing_uniqueness
	save(results_tandem_NX_annot_disc_missing_uniqueness_2UniqPep, file = paste0(projectDir, "Results/results_tandem_NX_annot_disc_missing_uniqueness_2UniqPep.rda"))
	results_tandem_NX_annot_disc_missing_uniqueness_f <- results_tandem_NX_annot_disc_missing_uniqueness_2UniqPep
	results_tandem_NX_annot_disc_missing_uniqueness_ff <- results_tandem_NX_annot_disc_missing_uniqueness_2UniqPep
}

#venn diagrams

pdf(file = paste0(projectDir, "Plots/results_Venn_MissingNXProteins_with_2_peptides_PerSearchEngine_after_uniqueness_checker.pdf"), width = 10, height = 10)
compare4List(paste(results_comet_NX_annot_disc_missing_uniqueness_2UniqPep$NextprotID), paste(results_mascot_NX_annot_disc_missing_uniqueness_2UniqPep$NextprotID), paste(results_omssa_NX_annot_disc_missing_uniqueness_2UniqPep$NextprotID), paste(results_tandem_NX_annot_disc_missing_uniqueness_2UniqPep$NextprotID), "Comet results", "Mascot results", "Omssa results", "Tandem results",  "ProtFDRGlobal - Missing proteins (NX) for each search engine with 2 peptides")
#compare4List(paste(results_comet_NX_annot_disc_missing_uniqueness_2UniqPep$ProteinID), paste(results_mascot_NX_annot_disc_missing_uniqueness_2UniqPep$ProteinID), paste(results_omssa_NX_annot_disc_missing_uniqueness_2UniqPep$ProteinID), paste(results_tandem_NX_annot_disc_missing_uniqueness_2UniqPep$ProteinID), "Comet results", "Mascot results", "Omssa results", "Tandem results",  "ProtFDRGlobal - Missing proteins (NX) for each search engine")
dev.off()

library(UpSetR)

listInput <- list(mascot = paste(results_mascot_NX_annot_disc_missing_uniqueness_2UniqPep$NextprotID) , comet = paste(results_comet_NX_annot_disc_missing_uniqueness_2UniqPep$NextprotID), omssa =paste(results_omssa_NX_annot_disc_missing_uniqueness_2UniqPep$NextprotID), tandem=paste(results_tandem_NX_annot_disc_missing_uniqueness_2UniqPep$NextprotID)  )


pdf(file = paste0(projectDir, "Plots/results_upsetR_MissingNXProteins_with_2_peptides_PerSearchEngine_after_uniqueness_checker.pdf"), width = 10, height = 10)
upset(fromList(listInput), order.by = "freq")
dev.off()




pdf(file = paste0(projectDir, "Plots/results_Venn_PeptidesOfMissingNXProteinsPerSearchEngine_with_2_peptides_after_uniqueness_checker.pdf"), width = 10, height = 10)
compare4List(paste(results_comet_NX_annot_disc_missing_uniqueness_2UniqPep$PeptideSeq), paste(results_mascot_NX_annot_disc_missing_uniqueness_2UniqPep$PeptideSeq), paste(results_omssa_NX_annot_disc_missing_uniqueness_2UniqPep$PeptideSeq), paste(results_tandem_NX_annot_disc_missing_uniqueness_2UniqPep$PeptideSeq), "Comet results", "Mascot results", "Omssa results", "Tandem results",  "ProtFDRGlobal - Peptides of missing proteins (NX) for each search engine with 2 peptides")
dev.off()

listInput <- list(mascot = paste(results_mascot_NX_annot_disc_missing_uniqueness_2UniqPep$PeptideSeq) , comet = paste(results_comet_NX_annot_disc_missing_uniqueness_2UniqPep$PeptideSeq), omssa =paste(results_omssa_NX_annot_disc_missing_uniqueness_2UniqPep$PeptideSeq), tandem=paste(results_tandem_NX_annot_disc_missing_uniqueness_2UniqPep$PeptideSeq)  )


pdf(file = paste0(projectDir, "Plots/results_upsetR_PeptidesOfMissingNXProteinsPerSearchEngine_with_2_peptides_after_uniqueness_checker.pdf"), width = 10, height = 10)
upset(fromList(listInput), order.by = "freq")
dev.off()


pdf(file=paste0(projectDir,"Plots/results_Venn_MissingNXProteins_Per_SearchEngine_one_hit_wonders_after_uniqueness_checker.pdf"), width = 10, height = 10)
compare4List(paste(results_comet_NX_annot_disc_missing_uniqueness_one_hit$NextprotID), paste(results_mascot_NX_annot_disc_missing_uniqueness_one_hit$NextprotID), paste(results_omssa_NX_annot_disc_missing_uniqueness_one_hit$NextprotID), paste(results_tandem_NX_annot_disc_missing_uniqueness_one_hit$NextprotID), "Comet results", "Mascot results", "Omssa results", "Tandem results",  "ProtFDRGlobal - Missing proteins (NX)-One hit wonders")
dev.off()

listInput <- list(mascot = paste(results_mascot_NX_annot_disc_missing_uniqueness_one_hit$NextprotID) , comet = paste(results_comet_NX_annot_disc_missing_uniqueness_one_hit$NextprotID), omssa =paste(results_omssa_NX_annot_disc_missing_uniqueness_one_hit$NextprotID), tandem=paste(results_tandem_NX_annot_disc_missing_uniqueness_one_hit$NextprotID)  )


pdf(file = paste0(projectDir, "Plots/results_upsetR_MissingNXProteins_Per_SearchEngine_one_hit_wonders_after_uniqueness_checker.pdf"), width = 10, height = 10)
upset(fromList(listInput), order.by = "freq")
dev.off()


pdf(file = paste0(projectDir, "Plots/results_Venn_PeptidesOfMissingNXProteinsPerSearchEngine_one_hit_wonders_after_uniqueness_checker.pdf"), width = 10, height = 10)
compare4List(paste(results_comet_NX_annot_disc_missing_uniqueness_one_hit$PeptideSeq), paste(results_mascot_NX_annot_disc_missing_uniqueness_one_hit$PeptideSeq), paste(results_omssa_NX_annot_disc_missing_uniqueness_one_hit$PeptideSeq), paste(results_tandem_NX_annot_disc_missing_uniqueness_one_hit$PeptideSeq), "Comet results", "Mascot results", "Omssa results", "Tandem results",  "ProtFDRGlobal - Peptides of missing proteins (NX) for each search engine with 1 peptides")
dev.off()

listInput <- list(mascot = paste(results_mascot_NX_annot_disc_missing_uniqueness_one_hit$PeptideSeq) , comet = paste(results_comet_NX_annot_disc_missing_uniqueness_one_hit$PeptideSeq), omssa =paste(results_omssa_NX_annot_disc_missing_uniqueness_one_hit$PeptideSeq), tandem=paste(results_tandem_NX_annot_disc_missing_uniqueness_one_hit$PeptideSeq)  )


pdf(file = paste0(projectDir, "Plots/results_upsetR_PeptidesOfMissingNXProteinsPerSearchEngine_one_hit_wonders_after_uniqueness_checker.pdf"), width = 10, height = 10)
upset(fromList(listInput), order.by = "freq")
dev.off()



# print data to table
# results_NX<-results_tandem_NX_annot_disc_missing_uniqueness_f
# results_NX_2<-results_tandem_NX_annot_disc_missing_uniqueness_ff
  results_NX <- rbind(results_comet_NX_annot_disc_missing_uniqueness_f, results_mascot_NX_annot_disc_missing_uniqueness_f, results_omssa_NX_annot_disc_missing_uniqueness_f, results_tandem_NX_annot_disc_missing_uniqueness_f)
  results_NX_2<-rbind(results_comet_NX_annot_disc_missing_uniqueness_ff, results_mascot_NX_annot_disc_missing_uniqueness_ff, results_omssa_NX_annot_disc_missing_uniqueness_ff, results_tandem_NX_annot_disc_missing_uniqueness_ff)
write.table(results_NX, file = paste0(projectDir, "Results/results_NX.txt"), col.names=TRUE, row.names = FALSE, sep="\t", quote=FALSE)
write.table(results_NX_2, file = paste0(projectDir, "Results/results_NX_FF.txt"), col.names=TRUE, row.names = FALSE, sep="\t", quote=FALSE)


##### MATRIZ DE MISSINGS CON 2 PEPTIDOS UNICOS

mascot_selected<-results_mascot_NX_annot_disc_missing[which(results_mascot_NX_annot_disc_missing$PSM %in% results_mascot_NX_annot_disc_missing_uniqueness_2UniqPep$PSM ),]
comet_selected<-results_comet_NX_annot_disc_missing[which(results_comet_NX_annot_disc_missing$PSM %in% results_comet_NX_annot_disc_missing_uniqueness_2UniqPep$PSM ),]
omssa_selected<-results_omssa_NX_annot_disc_missing[which(results_omssa_NX_annot_disc_missing$PSM %in% results_omssa_NX_annot_disc_missing_uniqueness_2UniqPep$PSM ),]
tandem_selected<-results_tandem_NX_annot_disc_missing[which(results_tandem_NX_annot_disc_missing$PSM %in% results_tandem_NX_annot_disc_missing_uniqueness_2UniqPep$PSM ),]

# quedarme con mejor hit para cada query
find_best_psm_score<-function(mascot_selected){
	mascot_best_psm<-list()
	mascot_best_score<-list()
	mascot_best_protein<-list()
	for (i in 1:length(unique(paste(mascot_selected$Query)))){
		query<-unique(paste(mascot_selected$Query))[i]
	    tmp<-mascot_selected[which(mascot_selected$Query== query),]
	    if((nrow(tmp) >=2 )& length(unique(paste(tmp$score) >1))){
	        tmp_selected<-tmp[which.max(as.numeric(paste(tmp$score))),]
	        mascot_best_psm[[i]]<-paste(tmp_selected$PSM)
			mascot_best_score[[i]]<-paste(tmp_selected$score)
			mascot_best_protein[[i]]<-paste(tmp_selected$NextprotID)
		}else if(nrow(tmp)==1){
	        mascot_best_psm[[i]]<-paste(tmp$PSM)
			mascot_best_score[[i]]<-paste(tmp$score)
			mascot_best_protein[[i]]<-paste(tmp$NextprotID)
	    }

	}

	mascot_best_psm_score<-unlist(mascot_best_score)
	mascot_best_psm<-unlist(mascot_best_psm)

	mascot_df<-data.frame("PSM"=paste(mascot_best_psm), "score_mascot"=paste(mascot_best_psm_score), "protein_mascot"=paste(mascot_best_protein))
	return(mascot_df)
}


mascot_df<-find_best_psm_score(mascot_selected)
comet_df<-find_best_psm_score(comet_selected)
names(comet_df)[2]<-"score_comet"
names(comet_df)[3]<-"protein_comet" #no hay
omssa_df<-find_best_psm_score(omssa_selected)
names(omssa_df)[2]<-"score_omssa"
names(omssa_df)[3]<-"protein_omssa"
tandem_df<-find_best_psm_score(tandem_selected)
names(tandem_df)[2]<-"score_tandem"
names(tandem_df)[3]<-"protein_tandem"

list_proteins<-c(paste(tandem_df$protein_tandem),paste(omssa_df$protein_omssa), paste(comet_df$protein_comet),paste(mascot_df$protein_mascot))
list_psm<-c(paste(tandem_df$PSM),paste(omssa_df$PSM), paste(comet_df$PSM),paste(mascot_df$PSM))

df_tmp<-data.frame("PSM"=list_psm,"Protein"=list_proteins)


mascot_df<-mascot_df[,c("PSM","score_mascot")]
comet_df<-comet_df[,c("PSM","score_comet")]
omssa_df<-omssa_df[,c("PSM","score_omssa")]
tandem_df<-tandem_df[,c("PSM","score_tandem")]





list.psm<-c(paste(tandem_df$PSM), paste(mascot_df$PSM), paste(omssa_df$PSM), paste(comet_df$PSM))
Query_df<-data.frame("PSM"=list.psm)

Query_df<-merge(Query_df, mascot_df, by.x="PSM", all.x=T)
Query_df<-unique(Query_df)
Query_df<-merge(Query_df, comet_df,by="PSM", all.x=T)
Query_df<-unique(Query_df)
Query_df<-merge(Query_df, tandem_df,by="PSM", all.x=T)
Query_df<-unique(Query_df)
Query_df<-merge(Query_df, omssa_df,by="PSM", all.x=T)
Query_df<-unique(Query_df)
Query_df<-merge(Query_df, df_tmp,by="PSM", all.x=T)
Query_df<-Query_df[order(Query_df$Protein),]

write.table(Query_df, file=paste0(projectDir,"Results/matrix_best_PSM_proteins_with_2Uniq_Peptides.txt"), col.names=T,row.names=F,quote=F,sep="\t")

mascot_selected<-results_mascot_NX_annot_disc_missing[which(results_mascot_NX_annot_disc_missing$PSM %in% results_mascot_NX_annot_disc_missing_uniqueness_one_hit$PSM ),]
comet_selected<-results_comet_NX_annot_disc_missing[which(results_comet_NX_annot_disc_missing$PSM %in% results_comet_NX_annot_disc_missing_uniqueness_one_hit$PSM ),]
omssa_selected<-results_omssa_NX_annot_disc_missing[which(results_omssa_NX_annot_disc_missing$PSM %in% results_omssa_NX_annot_disc_missing_uniqueness_one_hit$PSM ),]
tandem_selected<-results_tandem_NX_annot_disc_missing[which(results_tandem_NX_annot_disc_missing$PSM %in% results_tandem_NX_annot_disc_missing_uniqueness_one_hit$PSM ),]

# quedarme con mejor hit para cada query
mascot_df<-find_best_psm_score(mascot_selected)
comet_df<-find_best_psm_score(comet_selected)
names(comet_df)[2]<-"score_comet"
names(comet_df)[3]<-"protein_comet" #no hay
omssa_df<-find_best_psm_score(omssa_selected)
names(omssa_df)[2]<-"score_omssa"
names(omssa_df)[3]<-"protein_omssa"
tandem_df<-find_best_psm_score(tandem_selected)
names(tandem_df)[2]<-"score_tandem"
names(tandem_df)[3]<-"protein_tandem"

list_proteins<-c(paste(tandem_df$protein_tandem),paste(omssa_df$protein_omssa), paste(comet_df$protein_comet),paste(mascot_df$protein_mascot))
list_psm<-c(paste(tandem_df$PSM),paste(omssa_df$PSM), paste(comet_df$PSM),paste(mascot_df$PSM))

df_tmp<-data.frame("PSM"=list_psm,"Protein"=list_proteins)


mascot_df<-mascot_df[,c("PSM","score_mascot")]
comet_df<-comet_df[,c("PSM","score_comet")]
omssa_df<-omssa_df[,c("PSM","score_omssa")]
tandem_df<-tandem_df[,c("PSM","score_tandem")]





list.psm<-c(paste(tandem_df$PSM), paste(mascot_df$PSM), paste(omssa_df$PSM), paste(comet_df$PSM))
Query_df<-data.frame("PSM"=list.psm)

Query_df<-merge(Query_df, mascot_df, by.x="PSM", all.x=T)
Query_df<-unique(Query_df)
Query_df<-merge(Query_df, comet_df,by="PSM", all.x=T)
Query_df<-unique(Query_df)
Query_df<-merge(Query_df, tandem_df,by="PSM", all.x=T)
Query_df<-unique(Query_df)
Query_df<-merge(Query_df, omssa_df,by="PSM", all.x=T)
Query_df<-unique(Query_df)
Query_df<-merge(Query_df, df_tmp,by="PSM", all.x=T)
Query_df<-Query_df[order(Query_df$Protein),]

Query_df<-unique(Query_df)

write.table(Query_df, file=paste0(projectDir,"Results/matrix_best_PSM_proteins_with_one_hits.txt"), col.names=T,row.names=F,quote=F,sep="\t")


##### todos juntos

mascot_df<-find_best_psm_score(results_mascot_NX_annot_disc_missing)
comet_df<-find_best_psm_score(results_comet_NX_annot_disc_missing)
names(comet_df)[2]<-"score_comet"
names(comet_df)[3]<-"protein_comet" #no hay
omssa_df<-find_best_psm_score(results_omssa_NX_annot_disc_missing)
names(omssa_df)[2]<-"score_omssa"
names(omssa_df)[3]<-"protein_omssa"
tandem_df<-find_best_psm_score(results_tandem_NX_annot_disc_missing)
names(tandem_df)[2]<-"score_tandem"
names(tandem_df)[3]<-"protein_tandem"

list_proteins<-c(paste(tandem_df$protein_tandem),paste(omssa_df$protein_omssa), paste(comet_df$protein_comet),paste(mascot_df$protein_mascot))
list_psm<-c(paste(tandem_df$PSM),paste(omssa_df$PSM), paste(comet_df$PSM),paste(mascot_df$PSM))

df_tmp<-data.frame("PSM"=list_psm,"Protein"=list_proteins)


mascot_df_tmp<-mascot_df[,c("PSM","score_mascot")]
comet_df_tmp<-comet_df[,c("PSM","score_comet")]
omssa_df_tmp<-omssa_df[,c("PSM","score_omssa")]
tandem_df_tmp<-tandem_df[,c("PSM","score_tandem")]





list.psm<-c(paste(tandem_df$PSM), paste(mascot_df$PSM), paste(omssa_df$PSM), paste(comet_df$PSM))
Query_df<-data.frame("PSM"=list.psm)

Query_df<-merge(Query_df, mascot_df_tmp, by.x="PSM", all.x=T)
Query_df<-unique(Query_df)
Query_df<-merge(Query_df, comet_df_tmp,by="PSM", all.x=T)
Query_df<-unique(Query_df)
Query_df<-merge(Query_df, tandem_df_tmp,by="PSM", all.x=T)
Query_df<-unique(Query_df)
Query_df<-merge(Query_df, omssa_df_tmp,by="PSM", all.x=T)
Query_df<-unique(Query_df)
Query_df<-merge(Query_df, df_tmp,by="PSM", all.x=T)
Query_df<-Query_df[order(Query_df$Protein),]

Query_df<-unique(Query_df)

write.table(Query_df, file=paste0(projectDir,"Results/matrix_best_PSM_ALL.txt"), col.names=T,row.names=F,quote=F,sep="\t")

### Para seleccionarlos yo mejor los espectros

find_best_score<-function(dataframe_selected){
	dataframe_best_psm<-list()
	dataframe_best_score<-list()
	dataframe_best_protein<-list()
	for (i in 1:length(unique(paste(dataframe_selected$PeptideSeq)))){
		peptide<-unique(paste(dataframe_selected$PeptideSeq))[i]
	    tmp<-dataframe_selected[which(dataframe_selected$PeptideSeq== peptide),]
	    if((nrow(tmp) >=2 )& length(unique(paste(tmp$score) >1))){
	        tmp_selected<-tmp[which.max(as.numeric(paste(tmp$score))),]
	        dataframe_best_psm[[i]]<-paste(tmp_selected$PSM)
			dataframe_best_score[[i]]<-paste(tmp_selected$score)
			dataframe_best_protein[[i]]<-paste(tmp_selected$NextprotID)
		}else if(nrow(tmp)==1){
	        dataframe_best_psm[[i]]<-paste(tmp$PSM)
			dataframe_best_score[[i]]<-paste(tmp$score)
			dataframe_best_protein[[i]]<-paste(tmp$NextprotID)
	    }

	}

	dataframe_best_psm_score<-unlist(dataframe_best_score)
	dataframe_best_psm<-unlist(dataframe_best_psm)

	dataframe_df<-data.frame("PSM"=paste(dataframe_best_psm), "score_dataframe"=paste(dataframe_best_psm_score), "protein_dataframe"=paste(dataframe_best_protein))
	return(dataframe_df)
}





#####A PARTIR DE AQUI LO Q YO HE HECHO

# output_results<-results_NX_2[(results_NX_2$NofDiscPeptides>0),]
# write.table(output_results, file = paste0(projectDir, "Results/results_NX_filtered_NofDisc.txt"), col.names=TRUE, row.names = FALSE, sep="\t", quote=FALSE)
#
# output_results_for_agg<-output_results2[,c("PeptideSeq","NexprotID","Chr")]
#
# output_results_for_agg_unique<-unique(output_results_for_agg)
#
# load("/home/nostromo/data/pepe/02_neXtprot_20180117_Feb18/nextProtXChrXENSP.RData")
#
# db_selected<-nextProtXChrXENSP[,c("NextprotID","GeneName")]
#
# output_final<-merge(output_results_for_agg, db_selected, by.x = c("NexprotID"), by.y = c("NextprotID"), all.x = TRUE)
# output_final<-unique(output_final)
# write.table(output_final,file=paste0(projectDir, "Results/peptides_NX_all_SE.txt"), col.names=TRUE, row.names = FALSE, sep="\t", quote=FALSE))
