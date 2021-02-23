args=(commandArgs(TRUE))

xmlFile_in <- args[1]

cat(system("date"))
cat(xmlFile_in,"\n")

require(XML, quietly=T)
require(stringr,quietly=T)
require(gdata,quietly=T)

data<-xmlParse(xmlFile_in)
xml_data_raw <- xmlToList(data)
#save(xml_data_raw, file = paste(gsub(".xml","",xmlFile), "_raw.rda", sep = ""))
xml_data_unlist <- unlist(xml_data_raw)

parsedataXML<-function(xml_data_unlist){

    
    indexes <- which(names(xml_data_unlist) == "entry-list.entry.overview.protein-existence.value")
    cat("proteins =", length(indexes)-1,"\n")

    xml_data_df <- data.frame("NX_code"=as.character(),
    "GO_code"=as.character()
    )
    tests<-c()
    names<-c()

        for (i in 1:length(indexes)) {
        	if (i == 1) {
          # cat(i, "of ",length(indexes)-1,"\n")    
			index <- 1
			index_fin <- indexes[i] - 1
            cat(index,"/",index_fin,"\n")
		}else{
          #  cat(i, "of ",length(indexes)-1,"\n")
			index <- indexes[i-1]
			index_fin <- indexes[i] - 1
            cat(index,"/",index_fin,"\n")
		}
        tmp <- xml_data_unlist[index:index_fin]
        tmp_df <- data.frame("name" = names(tmp), "value" = tmp)
            if(length(unique(paste(tmp_df[tmp_df$name =="entry-list.entry..attrs.accession","value"]))) ==1 & (length(unique(paste(tmp_df[tmp_df$name == "entry-list.entry.annotation-list.annotation-category.annotation.cv-term..attrs.accession","value"]))) >=1)){
                
                subset_GO<-paste(tmp_df[tmp_df$name == "entry-list.entry.annotation-list.annotation-category.annotation.cv-term..attrs.accession","value"])
                subset_GO<-subset_GO[startsWith(subset_GO, "GO:")]
                subset_GO<-unique(subset_GO)
                if(length(subset_GO)==0){
                    subset_GO<-"NA"
                }

                df<-data.frame("NX_code"= rep(unique(paste(tmp_df[tmp_df$name =="entry-list.entry..attrs.accession","value"])), length(subset_GO)), "GO_code"=paste(subset_GO) )

                xml_data_df<-rbind(xml_data_df,df)

            }else if(length(unique(paste(tmp_df[tmp_df$name =="entry-list.entry..attrs.accession","value"]))) ==1 & (length(unique(paste(tmp_df[tmp_df$name == "entry-list.entry.annotation-list.annotation-category.annotation.cv-term..attrs.accession","value"]))) ==0))
            {
                
                
                subset_GO<-"NA"

                df<-data.frame("NX_code"= rep(unique(paste(tmp_df[tmp_df$name =="entry-list.entry..attrs.accession","value"])), length(subset_GO)), "GO_code"=paste(subset_GO) )

                xml_data_df<-rbind(xml_data_df,df)

            
            }

        }
        #### ultimo caso
        tmp<-xml_data_unlist[index_fin:length(xml_data_unlist)]
        tmp_df <- data.frame("name" = names(tmp), "value" = tmp)

            if(length(unique(paste(tmp_df[tmp_df$name =="entry-list.entry..attrs.accession","value"]))) ==1 & (length(unique(paste(tmp_df[tmp_df$name == "entry-list.entry.annotation-list.annotation-category.annotation.cv-term..attrs.accession","value"]))) >=1))
            {

                subset_GO<-paste(tmp_df[tmp_df$name == "entry-list.entry.annotation-list.annotation-category.annotation.cv-term..attrs.accession","value"])
                subset_GO<-subset_GO[startsWith(subset_GO, "GO:")]
                subset_GO<-unique(subset_GO)
                if(length(subset_GO)==0){
                    subset_GO<-"NA"
                }

                df<-data.frame("NX_code"= rep(unique(paste(tmp_df[tmp_df$name =="entry-list.entry..attrs.accession","value"])), length(subset_GO)), "GO_code"=paste(subset_GO) )

                xml_data_df<-rbind(xml_data_df,df)

            }else if(length(unique(paste(tmp_df[tmp_df$name =="entry-list.entry..attrs.accession","value"]))) ==1 & (length(unique(paste(tmp_df[tmp_df$name == "entry-list.entry.annotation-list.annotation-category.annotation.cv-term..attrs.accession","value"]))) ==0))
            {
                subset_GO<-"NA"

                df<-data.frame("NX_code"= rep(unique(paste(tmp_df[tmp_df$name =="entry-list.entry..attrs.accession","value"])), length(subset_GO)), "GO_code"=paste(subset_GO) )

                xml_data_df<-rbind(xml_data_df,df)

            }
        

    return((xml_data_df)

}

Nextprot_output<-parsedataXML(xml_data_unlist)

dir_our<-dirname(xmlFile_in)
name_s<-basename(xmlFile_in)
name_out<-lapply(strsplit(paste(name_s),"\\."),"[",1)

write.table(Nextprot_output, file=paste0(dir_our,"/",name_out,".tsv"), col.names=T, row.names=F, quote=F, sep="\t")

cat(system("date"),"\n")

cat("DONE", name_s, "\n")



