#!/usr/bin/Rscript

args=(commandArgs(TRUE))

OUTPUT_DIR <- args[1] #file Shotgun___.txt
method<-args[2] #dot, dist, cos, mutual
minn_AA<-args[3] # minimum peptideseq

library(stringr)
library(Biostrings)
library(doBy)
source("/home/margaret/data/01_Rscripts/A_Funciones/funcionesShotgun.R")
cat("Starting the script \n")

psmFDRvalue = 0.01
protFDRvalue = 0.01
pepScore="score"
#decoy_id="decoy"
databaseLoaded <- 0
maxn_AA<-50
concat_decoy = 0

OUTPUT_DIR_TARGET<-paste0(OUTPUT_DIR,"TARGET/")
OUTPUT_DIR_DECOY<-paste0(OUTPUT_DIR,"DECOY/")

filestarget<-list.files(OUTPUT_DIR_TARGET, pattern=".Rdata")
filestarget<-paste0(OUTPUT_DIR_TARGET,filestarget)

filesdecoy<-list.files(OUTPUT_DIR_DECOY, pattern=".Rdata")
filesdecoy<-paste0(OUTPUT_DIR_DECOY,filesdecoy)

dataPSMMat<-data.frame()

for(i in 1:length(filestarget)){
	loaded<-get(load(filestarget[i]))
	dataPSMMat<-rbind(dataPSMMat,loaded)
}


for(i in 1:length(filesdecoy)){
	loaded<-get(load(filesdecoy[i]))
	dataPSMMat<-rbind(dataPSMMat,loaded)
}

dataPSMMat$PeptideSeq<-lapply(strsplit(paste(dataPSMMat$Name),"/"),"[",1)
dataPSMMat$PeptideSeq<-paste(dataPSMMat$PeptideSeq)

cargar

readMSPFile<-function(x){
    con <- file(x, open="r")
    Names_spectra_ALL<-list()
    MW_ALL<-list()
    Comments_ALL<-list()
    NumPeaks_ALL<-list()
    Peaks_ALL<-list()
    PrecursorMZ_ALL<-list()
    counter<-0
    line <- 'Hello'
    reading <- T
    while (reading == T){
        line <- readLines(con, n=1)
        #print(reading)
        if(length(line) == 0){break}
        if(grepl("^Name", line)){
            counter = counter+1
            Names_spectra_ALL[counter]<-line
        }else if(grepl("^MW",line)){
            MW_ALL[counter]<-line
        }else if(grepl("^PrecursorMZ",line)){
            PrecursorMZ_ALL[counter]<-line
        }else if(grepl("^Comment",line)){
            Comments_ALL[counter]<-line
        }else if(grepl("^Num", line)){
            NumPeaks_ALL[counter]<-line
        }else if(grepl('^[0-9]', line)){
            if(is.null(Peaks_ALL[counter][1])){
                Peaks_ALL[counter]<-line
            }else{
                Peaks_ALL[counter]<-paste(Peaks_ALL[counter], '##',line)
            }
         }
    }
    close(con)
    Names_spectra_ALL<-sapply(strsplit(paste(Names_spectra_ALL),": "),"[",2)
    NumPeaks_ALL<-as.numeric(str_extract(NumPeaks_ALL,'[0-9\\.]+'))
    MW_All<-as.numeric(str_extract(MW_ALL,'[0-9\\.]+'))
    #Charge_all<-sapply(strsplit(paste(Names_spectra_ALL),"\\"),"[",2)
    #Charge_all<-as.numeric(str_extract(Charge_all,'[0-9\\.]+' ))
    PrecursorMZ_ALL<-as.numeric(str_extract(PrecursorMZ_ALL,'[0-9\\.]+'))
    Peaks_ALL<-sapply(strsplit(paste(Peaks_ALL),"NULL ##"),"[",2)

    ProteinNameTMP<-sapply(strsplit(paste(Comments_ALL),"NISTProtein="),"[",2)
    ProteinName<-sapply(strsplit(paste(ProteinNameTMP),"\\("),"[",1)
	#ProteinName<-sapply(strsplit(paste(ProteinNameTMP),"\\("),"[",1)

	#limpiar campos y generar df


    MSP_df <-data.frame("Name" = paste(Names_spectra_ALL), "Lib_ID" = paste(seq(1:length(Names_spectra_ALL))), "MW" = paste(MW_All), "PrecursorMZ" = paste(PrecursorMZ_ALL), "NumPeaks" = paste(NumPeaks_ALL), "Comments"=paste(Comments_ALL), "Peaks"=paste(Peaks_ALL), "ProteinName"=paste(ProteinName))
    return(MSP_df)

}






dataPSMMatXPROT2<-merge(dataPSMMatXPROT,all_peptides_TMP2, by.x="PeptideSeq", by.y="peptide", all.x=T)

cat("Filtering AA length ...\n")
thr_AA <- as.numeric(paste(minn_AA)) - 1
dataPSMMat_filterAA <- dataPSMMat[nchar(paste(dataPSMMat$PeptideSeq)) > thr_AA,]
thr_AA <- as.numeric(paste(maxn_AA)) + 1
dataPSMMat_filterAA <- dataPSMMat_filterAA[nchar(paste(dataPSMMat_filterAA$PeptideSeq)) < thr_AA,]

save(dataPSMMat_filterAA, file=paste(currentDir, "/results_Peptides_", search_engine_str, "_", dataset, "_dataPSMMat_filterAA.rda", sep = ""))
cat("Calculating psmFDR and filtering ...\n")
dataPSMMat_filterAA_psmFDR <- psmFDR(dataPSMMat_filterAA, pepScore=pepScore, decoy_id=decoy_id, concat_decoy = concat_decoy, prot_id = prot_id)
dataPSMMat_filterAA_psmFDR_Filter <- dataPSMMat_filterAA_psmFDR[dataPSMMat_filterAA_psmFDR$psmFDR < psmFDRvalue,]
save(dataPSMMat_filterAA_psmFDR_Filter,file=paste(currentDir, "/results_Peptides_", search_engine_str, "_", dataset, "_dataPSMMat_filterAA_PSMFDR_filter.rda", sep = "") )


cat("Calculating pepFDR and filtering ...\n")
dataPSMMat_filterAA_psmFDR_Filter_pepFDR <- pepFDR(dataPSMMat_filterAA_psmFDR_Filter, pepScore=pepScore, concat_decoy = concat_decoy, pep_id = pep_col)
dataPSMMat_filterAA_psmFDR_Filter_pepFDR_Filter <- dataPSMMat_filterAA_psmFDR_Filter_pepFDR[dataPSMMat_filterAA_psmFDR_Filter_pepFDR$pepFDR < 0.01,]
save(dataPSMMat_filterAA_psmFDR_Filter_pepFDR_Filter, file = paste(currentDir, "/results_Peptides_", search_engine_str, "_", dataset, "_dataPSMMat_filterAA_psmFDR_Filter_pepFDR_Filter.rda", sep = ""))
cat("Calculating protFDR and filtering ...\n")
dataPSMMat_filterAA_psmFDR_Filter_pepFDR_Filter_protFDR <- protFDR(dataPSMMat_filterAA_psmFDR_Filter_pepFDR_Filter, pepScore=pepScore, concat_decoy = concat_decoy, prot_id = prot_id)
dataPSMMat_filterAA_psmFDR_Filter_pepFDR_Filter_protFDR_Filter <- dataPSMMat_filterAA_psmFDR_Filter_pepFDR_Filter_protFDR[dataPSMMat_filterAA_psmFDR_Filter_pepFDR_Filter_protFDR$protFDR < 0.01,]
save(dataPSMMat_filterAA_psmFDR_Filter_pepFDR_Filter_protFDR_Filter, file=paste(currentDir, "/results_Peptides_", search_engine_str, "_", dataset, "_dataPSMMat_filterAA_psmFDR_Filter_pepFDR_Filter_protFDR_Filter.rda", sep = ""))
