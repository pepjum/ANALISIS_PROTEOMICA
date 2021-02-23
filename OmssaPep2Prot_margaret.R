####################################################################################
###
args=(commandArgs(TRUE))


### Rscript TandemPep2Prot_mlsanchez.R /home/nostromo/data/mlsanchez/proteomica/ejemplo_experimento/ Exp1_M3 /home/nostromo/data/mlsanchez/proteomica/ejemplo_experimento/Tandem_Files/Shotgun_tandem_Exp1_M3.txt uniprot_sprot_2017_12_CRAP.fasta DECOY 0


# params
currentFolder <- args[1]
# currentFolder <- "/mnt/beegfs/agarin/dato-activo/03_Analysis/agarin/26_Navarrabiomed_Missing_Ene18/"
dataset <- args[2]
# dataset <- "349-028"
 txtFileName <- args[3]
# txtFileName <- "/mnt/beegfs/agarin/dato-activo/03_Analysis/agarin/26_Navarrabiomed_Missing_Ene18/Tandem_Files/Shotgun_tandem_349-028.txt"
DBtargetFile <- args[4]
# DBtargetFile <- "uniprot_sprot_2017_12.fasta"
decoy_id <- args[5]
# decoy_id <- "DECOY"
todo <- args[6]
# TOD; 1 si solo hay que hacer target; 2 si solo hay que hacer decoy; 0 si hay que hacer los dos.

library(data.table)
library(Biostrings)
library(doBy)
source("/home/margaret/data/01_Rscripts/A_Funciones/funcionesShotgun.R")

mgfFolder <- paste(currentFolder, "MGFFiles", sep = "")
tandemFolder <- paste(currentFolder, "Omssa_Files","", sep = "")
datasetMgfFolder <- paste(mgfFolder, dataset, sep = "/")
datasetTandemTargetFolder <- paste(tandemFolder, dataset, "", sep = "/")
datasetTandemDecoyFolder <- paste(tandemFolder,"/", dataset, "-D/",sep = "")
mgfFiles <- dir(path = datasetMgfFolder, pattern = ".mgf")
DBtargetFile <- paste(currentFolder, DBtargetFile, sep = "")
dbName <- strsplit(basename(DBtargetFile), "\\.")[[1]][1]

todo <- as.numeric(paste(todo))

if (todo == 0 || todo == 1) {
	targetTandemFiles <- dir(path = datasetTandemTargetFolder, pattern = ".txt")
	targetTandemFiles <- targetTandemFiles[!grepl("peptides.tsv", targetTandemFiles)]
	targetTandemFiles <- targetTandemFiles[!grepl("_peptides_out.tsv", targetTandemFiles)]
	targetTandemFiles <- targetTandemFiles[!grepl("_corrected.tsv", targetTandemFiles)]
	targetTandemFiles <- paste(datasetTandemTargetFolder, targetTandemFiles, sep = "")
	targetOutTandemFiles <- paste(substr(targetTandemFiles, 1, nchar(targetTandemFiles) - 4), "_peptides.tsv", sep="")
  targetPMOutTandemFiles <- paste(substr(targetTandemFiles, 1, nchar(targetTandemFiles) - 4), "_peptides_out.tsv", sep="")
	targetCorrectedTandemFiles <- paste(substr(targetTandemFiles, 1, nchar(targetTandemFiles) - 4), "_corrected.tsv", sep="")
  targetLogTandemFiles <-  paste(substr(targetTandemFiles, 1, nchar(targetTandemFiles) - 4), "_peptides_log.log", sep="")
  for(i in 1:length(targetTandemFiles)) {

		tandemResult <- read.table(targetTandemFiles[i], header = TRUE, sep = "\t")
		write.table(tandemResult$peptide, file = targetOutTandemFiles[i], quote = FALSE, row.names = FALSE, col.names = FALSE, sep="\t")

			y <- paste('java -XX:ParallelGCThreads=4 -jar /opt/PeptideMatchCMD/PeptideMatchCMD_1.0.jar -a query -i ', currentFolder, '/PeptideMatch_', dbName, '_index -Q ',targetOutTandemFiles[i],' -l -e -o ', targetPMOutTandemFiles[i], ' > ', targetLogTandemFiles[i], sep = "")
			system(y)
			y <- NULL

			cat(targetPMOutTandemFiles[i])
			peptidesMatched <- read.csv2(targetPMOutTandemFiles[i], header = FALSE, sep = "\t", skip = 2)
			peptidesMatched <- peptidesMatched[,1:2]
			peptidesMatched <- unique(peptidesMatched)
			colnames(peptidesMatched) <- c("peptide", "protein_id")
			tandemResultMatched <- merge(tandemResult, peptidesMatched, by.x = "peptide", by.y = "peptide", all.x = TRUE)
			tandemResultMatched <- unique(tandemResultMatched)
			tandemResultMatched$Protein <- NULL
			tandemResultMatched$Protein <- tandemResultMatched$protein_id
			tandemResultMatched$protein_id <- NULL

			write.table(tandemResultMatched, file = targetCorrectedTandemFiles[i], quote = FALSE, row.names = FALSE, col.names = TRUE, sep="\t")

			tandemResult <- NULL
			tandemResultMatched <- NULL
			peptidesMatched <- NULL

	}
}

if (todo == 0 || todo == 2) {
	decoyTandemFiles <- dir(path = datasetTandemDecoyFolder, pattern = ".txt")
	decoyTandemFiles <- decoyTandemFiles[!grepl("peptides.tsv", decoyTandemFiles)]
	decoyTandemFiles <- decoyTandemFiles[!grepl("_peptides_out.tsv", decoyTandemFiles)]
	decoyTandemFiles <- decoyTandemFiles[!grepl("_corrected.tsv", decoyTandemFiles)]
	decoyTandemFiles <- paste(datasetTandemDecoyFolder, decoyTandemFiles, sep = "")
	decoyOutTandemFiles <- paste(substr(decoyTandemFiles, 1, nchar(decoyTandemFiles) - 4), "_peptides.tsv", sep="")
	decoyPMOutTandemFiles <- paste(substr(decoyTandemFiles, 1, nchar(decoyTandemFiles) - 4), "_peptides_out.tsv", sep="")
  decoyCorrectedTandemFiles <- paste(substr(decoyTandemFiles, 1, nchar(decoyTandemFiles) - 4), "_corrected.tsv", sep="")
  decoyLogTandemFiles <-  paste(substr(decoyTandemFiles, 1, nchar(decoyTandemFiles) - 4), "_peptides_log.log", sep="")
  for(i in 1:length(decoyTandemFiles)) {
		tandemResult <- read.table(decoyTandemFiles[i], header = TRUE, sep = "\t")
		write.table(tandemResult$PeptideSeq, file = decoyOutTandemFiles[i], quote = FALSE, row.names = FALSE, col.names = FALSE, sep="\t")

		y <- paste('java -XX:ParallelGCThreads=4 -jar /opt/PeptideMatchCMD/PeptideMatchCMD_1.0.jar -a query -i ', currentFolder, '/PeptideMatch_', dbName, '_', decoy_id,'_index -Q ',decoyOutTandemFiles[i],' -l -e -o ',decoyPMOutTandemFiles[i], ' > ', decoyLogTandemFiles[i], sep = "")
		system(y)
		y <- NULL
  	    peptidesMatched <- read.csv2(decoyPMOutTandemFiles[i], header = FALSE, sep = "\t")
		peptidesMatched <- peptidesMatched[,1:2]
		peptidesMatched <- unique(peptidesMatched)
		colnames(peptidesMatched) <- c("peptide", "protein_id")
		tandemResultMatched <- merge(tandemResult, peptidesMatched, by.x = "peptide", by.y = "peptide", all.x = TRUE)
		tandemResultMatched <- unique(tandemResultMatched)

		tandemResultMatched$Protein <- NULL
		tandemResultMatched$Protein <- tandemResultMatched$protein_id
		tandemResultMatched$protein_id <- NULL

		write.table(tandemResultMatched,  file = decoyCorrectedTandemFiles[i], quote = FALSE, row.names = FALSE, col.names = TRUE, sep="\t")

		tandemResult <- NULL
		tandemResultMatched <- NULL
		peptidesMatched <- NULL

	}

}

targetTandemFiles <- dir(path = datasetTandemTargetFolder, pattern = "_corrected.tsv")
targetTandemFiles <- targetTandemFiles[order(targetTandemFiles, decreasing=FALSE)]
targetTandemFiles <- paste(datasetTandemTargetFolder, targetTandemFiles, sep = "")

decoyTandemFiles <- dir(path = datasetTandemDecoyFolder, pattern = "_corrected.tsv")
decoyTandemFiles <- decoyTandemFiles[order(decoyTandemFiles, decreasing=FALSE)]
decoyTandemFiles <- paste(datasetTandemDecoyFolder, decoyTandemFiles, sep = "")

tmp1 <- data.frame(a = targetTandemFiles, b = "T", c = dataset, d = 1:length(targetTandemFiles))
tmp2 <- data.frame(a = decoyTandemFiles, b = "D", c = dataset, d = 1:length(decoyTandemFiles))

tmp <- rbind(tmp1, tmp2)
write.table(tmp, file=txtFileName, col.names = FALSE, quote = FALSE, sep = " ", row.names = FALSE)
