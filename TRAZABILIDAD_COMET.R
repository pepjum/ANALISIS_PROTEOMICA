args=(commandArgs(TRUE))
#### TRAZABILIDAD. SE PUEDE HACER DESDE CUALQUIER rDATA RESULTANTE DEL ANALISIS SHOTGUN. NO ES NECESARIO INTERCALARLO

MGF_folder<-args[1]
#out_file<-args[2] #"/home/jgonzalez.69/dato-activo/03_Analysis/jgonzalez69/13_EMBRIO_NOV18/Comet_files/trazability.Rdata"

samples<-list.dirs(MGF_folder)
samples<-samples[-1]

library(gdata)
library(stringr)
#trazability_out<-data.frame()

for(i in 1:length(samples)){

    files<-list.files(samples[i], pattern=".mgf")
    files_out<-paste0(lapply(strsplit(files,"\\."),"[",1),".index")
    files<-paste(samples[i], files, sep="/")
    files_out<-paste(samples[i], files_out, sep="/")

    for (j in 1:length(files)){

        y<- paste0('grep "^TITLE" ', files[j],' > ', files_out[j])
        cat(y,"\n")
        system(y)
    }

    # for(k in 1:length(files_out)){
    #     y<-paste0('sed -i -e "s/ /_/g" ', files_out[k])
    #     cat(k,"\n")
    #     system(y)
    # }

    files_index_comet<-list.files(samples[i], pattern=".index")
    files_index_comet<-paste(samples[i],files_index_comet, sep="/")

    #extra_trazability<-data.frame()
    df_tmp_final<-data.frame()
    for(j in 1:length(files_index_comet)){
        cat(files_index_comet[j],"\n")
        trazability<-read.table(files_index_comet[j], header=F, comment.char ='$')
        #fraction_tmp<-paste(lapply(strsplit(paste(files_index_comet[j]),"\\/"),"[", 10))
        fraction_tmp<-basename(files_index_comet[j]) #revisa esto para cada experimento. Tiene q ser el nombre de la fraccion - la extension index
        fraction_tmp2<-substr(fraction_tmp,1,nchar(fraction_tmp)-6)  #normalmente tiene que ser -6

        # if(startsWith(paste(trazability$V1[1]), "TITLE=Locus")){
        #     kk<-sapply(strsplit(paste(trazability$V1),"\\."),"[",4)
        #     df_tmp<-data.frame("scanstring"=kk)
        # }else if(startsWith(paste(trazability$V1[1]), "TITLE=Cmpd")){
        #     kk<-str_extract(paste(trazability$V6),"[0-9]+")
        #     df_tmp<-data.frame("scanstring"=kk)
        # }else{
         df_tmp<-data.frame("scanstring"=trazability[,1])
        # }
        #df_tmp$scanstring<-lapply(strsplit(paste(df_tmp$scanstring),"\\."),"[",4)
        df_tmp$scanstring<-lapply(strsplit(paste(df_tmp$scanstring),"="),"[",2)
        names(df_tmp)<-"scanstring"
        df_tmp$mgfIndexes<-seq(1:nrow(df_tmp))
        df_tmp$fraction_t<-fraction_tmp2
        df_tmp$dummyCode<-paste(paste(as.character(df_tmp$fraction_t)),paste(as.character(df_tmp$mgfIndexes)),sep="_")
        df_tmp<-df_tmp[,c(1,4)]
        df_tmp_final<-rbind(df_tmp_final,df_tmp)
        #trazability_out<-rbind(trazability_out, df_tmp)
    }
    tmp3<-dirname(files_index_comet[j])
    experiment_name<-basename(tmp3)
    mgf_folder<-dirname(tmp3)
    Comet_folder<-gsub("MGFFiles", "Comet_files", mgf_folder)
    out_file<-paste0("trazability_",experiment_name,".Rdata")
    out_file<-paste(Comet_folder,out_file,sep="/")
    #fraction_tmp2<-sapply(strsplit(paste(fraction_tmp),"\\."),"[",1)
    cat("saving", out_file,"\n")
    save(df_tmp_final, file=out_file)
    #extra_trazability<-rbind(extra_trazability,df_tmp)
}
cat("end")
