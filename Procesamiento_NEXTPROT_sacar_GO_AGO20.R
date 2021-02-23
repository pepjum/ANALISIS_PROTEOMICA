files<-list.files(pattern="xml.gz")

for(file in files){
    cat(file, "\n")
    y<-paste0("gzip -d ", file )
    system(y)
}


files<-list.files("/home/nostromo/data/pepe/21_neXtprot_2020_01_17/NEXTPROT_XML/xml",pattern="xml")

a<-paste0(getwd(),"/")

files<-paste0(a,files)

for(file in files){
    cat(file, "\n")
    y<-paste0("Rscript ~/data/pepe/scripts_tereshkova/ANALISIS_PROTEOMICA/parseNEXTPROT_XML_AGO20.R ", file )
    cat(y,"\n")
    system(y)
}


files<-list.files(pattern=".tsv)

nextprot_all<-data.frame()
for(file in files){
    cat(file"\n")
    file_loaded<-read.table(file, header=T)
    nextprot_all<-rbind(nextprot_all, file_loaded)
}

