#### NEXTPROT uPEs


library(dplyr)

#nextprot_enspxensg<-get(load("/home/nostromo/data/pepe/21_neXtprot_2020_01_17/nextProtXChrXENSPXENSG.RData"))
nextprot_enspxensg<-get(load("/home/nostromo/data/pepe/22_neXtProt_2020_7_17/nextProtXChrXENSPXENSG.RData"))


nextprot_enspxensg_sel<-nextprot_enspxensg[-which(nextprot_enspxensg$PE=="PE5"),]


"parseENCODE_eli" <- function(x) {

    tmp <- unlist(strsplit(unlist(strsplit(x, ";")), " "))
    tmp <- tmp[tmp != ""]
    tmp <- tmp[c(grep("gene_id",tmp)+1,grep("gene_type",tmp)+1,grep("gene_name",tmp)+1,grep("level",tmp)+1, ifelse(sum(grep("havana_gene",tmp))>0, grep("havana_gene",tmp)+1, NA))]
    return(tmp)

}

genecodev28_all <- read.table("~/data/00_References/gencode_gtf/gencode.v28.annotation.gtf", skip = 5, header = FALSE, sep = "\t")
genecodev28 <- genecodev28_all[genecodev28_all$V3 == "gene", ]
colnames(genecodev28) <- c("chr","DB","Type", "start","end","","strand","","Description")
genecodev28_tmp <- apply(as.data.frame(genecodev28[,9]), 1, parseENCODE_eli)
genecodev28_tmp <- t(genecodev28_tmp)
colnames(genecodev28_tmp) <- c("gene_id", "gene_type", "gene_name", "level", "havana_gene")
genecodev28_Annot <- data.frame(genecodev28_tmp, genecodev28[,-c(6,8:9)])
rownames(genecodev28_Annot) <- genecodev28_Annot[,1]
genecodev28_Annot$ENSG <- substr(genecodev28_Annot[,1],1,15)


nextprotXensgXgeneCode <- merge(nextprot_enspxensg_sel, genecodev28_Annot, by.x=9, by.y=12, all.x=F, all.y=F)
nextprotXensgXgeneCode<-nextprotXensgXgeneCode[,c("NextprotID","ENSG","Chr","GeneName","ProteinEvidence","PE","Missing","gene_id","gene_type","start","end","strand")]

#unknown_function<-read.table("/home/nostromo/data/pepe/21_neXtprot_2020_01_17/uPE_nextprot_2020-1-17.txt")
unknown_function<-read.table("/home/nostromo/data/pepe/22_neXtProt_2020_7_17/uPE_nextprot_2020_07_17.txt")



nextprotXensgXgeneCode$unknown_function<-nextprotXensgXgeneCode$NextprotID %in% paste(unknown_function$V1)*1

nextprotXensgXgeneCode_u  <- unique(nextprotXensgXgeneCode) %>% group_by(ENSG) %>% summarise(NextprotID=paste(unique(NextprotID), collapse=";"), GeneName=paste(unique(GeneName), collapse=";"), gene_id=paste(unique(gene_id), collapse=";"), gene_type=paste(unique(gene_type), collapse=";"), Chr=paste(unique(Chr), collapse=";"), start=paste(unique(start), collapse=";"), end=paste(unique(end), collapse=";"), ProteinEvidence=paste(unique(ProteinEvidence), collapse=";"), PE=paste(unique(PE), collapse=";"), Missing=paste(unique(Missing), collapse=";"), strand=paste(unique(strand),collapse=";"), unknown_function=paste(unique(unknown_function), collapse=";")) %>% as.data.frame

#save(nextprotXensgXgeneCode_u, file="~/data/jgonzalez/margaret/21_neXtprot_2020_01_17/nextprot_all_genecode+unknown.Rdata")
save(nextprotXensgXgeneCode_u, file="/home/nostromo/data/pepe/22_neXtProt_2020_7_17/nextprot_all_genecode+unknown.Rdata")

#write.table(nextprotXensgXgeneCode_u, file="~/data/jgonzalez/margaret/21_neXtprot_2020_01_17/nextprot_all_genecode+unknown.txt", col.names=T, row.names=F, quote=F, sep="\t")
write.table(nextprotXensgXgeneCode_u, file="/home/nostromo/data/pepe/22_neXtProt_2020_7_17/nextprot_all_genecode+unknown.txt", col.names=T, row.names=F, quote=F, sep="\t")



all_PE1<-nextprotXensgXgeneCode_u[which(nextprotXensgXgeneCode_u$PE=="PE1"),]

#save(all_PE1, file="~/data/jgonzalez/margaret/21_neXtprot_2020_01_17/all_PE1.Rdata")
#write.table(all_PE1, file="~/data/jgonzalez/margaret/21_neXtprot_2020_01_17/all_PE1.txt", col.names=T, row.names=F, quote=F, sep="\t")
write.table(all_PE1, file="/home/nostromo/data/pepe/22_neXtProt_2020_7_17/all_PE1.txt", col.names=T, row.names=F, quote=F, sep="\t")




PE1_unknown_function<-all_PE1[which(all_PE1$unknown_function==1),]
#save(PE1_unknown_function, file="~/data/jgonzalez/margaret/21_neXtprot_2020_01_17/uPE1_function.Rdata")

#write.table(PE1_unknown_function, file="~/data/jgonzalez/margaret/21_neXtprot_2020_01_17/uPE1_function.txt", col.names=T, row.names=F, quote=F, sep="\t")
write.table(PE1_unknown_function, file="~/home/nostromo/data/pepe/22_neXtProt_2020_7_17/uPE1_function.txt", col.names=T, row.names=F, quote=F, sep="\t")


PE1_known_function<-all_PE1[which(all_PE1$unknown_function==0),]
#save(PE1_known_function, file="~/data/jgonzalez/margaret/21_neXtprot_2020_01_17/PE1_known_function.Rdata")
#write.table(PE1_known_function, file="~/data/jgonzalez/margaret/21_neXtprot_2020_01_17/PE1_known_function.txt", col.names=T, row.names=F, quote=F, sep="\t")
write.table(PE1_known_function, file="/home/nostromo/data/pepe/22_neXtProt_2020_7_17/PE1_known_function.txt", col.names=T, row.names=F, quote=F, sep="\t")


