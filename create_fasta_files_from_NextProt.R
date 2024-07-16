# create fasta file from md5 file of nextProt
args=(commandArgs(TRUE))
library(BioStrings)

source("/home/margaret/data/01_Rscripts/A_Funciones/funcionesShotgun.R")

file_name<-args[1]

data<-read.table(file_name, header=T)

fileout<-paste0(dirname(file_name),"/nextProt.fasta")
seq<-paste(as.character(data$sequence))
names(seq) <- paste0(data$isoform)
sp_fa <- AAStringSet(seq)
writeXStringSet(sp_fa, fileout)



seq <- paste(as.character(data$sequence))
	#create decoy database (independent, without being concatenated)
cat("Creating decoy database \n")
seq_decoy <- sapply(seq,  FUN=function(x) .strPseudoreverseDecoy(x))
fileout_decoy<-paste0(dirname(file_name),"/nextProt_DECOY.fasta")
#write.fasta(paste(seq_decoy), paste0("DECOY",data$isoform), fileout_decoy, open = "w", nbchar = 60, as.string = FALSE)

names(seq_decoy) <- paste0("DECOY_",data$isoform)
sp_decoy_fa <- AAStringSet(seq_decoy)
cat("Writing the decoy database \n")
writeXStringSet(sp_decoy_fa, fileout_decoy)
