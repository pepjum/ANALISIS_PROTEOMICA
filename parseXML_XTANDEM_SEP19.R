#!/usr/bin/Rscript

args=(commandArgs(TRUE))

xmlFile_in <- args[1]
tsv_out<-args[2]

cat(system("date"))
cat(xmlFile_in,"\n")

require(XML)
require(stringr)

data<-xmlParse(xmlFile_in)
xml_data_raw <- xmlToList(data)
#save(xml_data_raw, file = paste(gsub(".xml","",xmlFile), "_raw.rda", sep = ""))
xml_data_unlist <- unlist(xml_data_raw)

parsedataXML<-function(xml_data_unlist){

    indexes <- which(names(xml_data_unlist) == "group..attrs.act")

    xml_data_df <- data.frame("Protein"=as.character(),
    "PeptideSeq"=as.character(),
    "scan"=as.character(),
    "PrecursorMZ"=as.character(),
    "charge"=as.character(),
    "expect"=as.character(),
    "index_mgf"=as.character(),
    "modification"=as.character(),
    "start"=as.character(),
    "end"=as.character(),
    "hyperscore"=as.character(),
    "missed_cleavages"=as.character())

    for (i in 1:length(indexes)) {
        	if (i == 1) {
			index <- 1
			index_fin <- indexes[i] - 1
		}else{
			index <- indexes[i-1]
			index_fin <- indexes[i] - 1
		}

    tmp <- xml_data_unlist[index:index_fin]
	tmp_df <- data.frame("name" = names(tmp), "value" = tmp)


        if ((length(paste(tmp_df[tmp_df$name == "group.protein.note.text","value"])) == 1) & (length(unique(paste(tmp_df[tmp_df$name == "group.protein.peptide.domain..attrs.seq","value"])))==1) & (length(unique(paste(tmp_df[tmp_df$name == "group.protein.peptide.domain..attrs.start","value"])))==1) & (length(unique(paste(tmp_df[tmp_df$name == "group.protein.peptide.domain..attrs.end","value"]))) == 1)){

            	df<- data.frame("Protein"=paste(tmp_df[tmp_df$name == "group.protein.note.text","value"]),
    				 "PeptideSeq"=unique(paste(tmp_df[tmp_df$name == "group.protein.peptide.domain..attrs.seq","value"])),
    				 "scan"=paste(tmp_df[tmp_df$name == "group.group.note.text","value"]),
    				 "PrecursorMZ"=paste(tmp_df[tmp_df$name == "group..attrs.mh","value"]),
    				 "charge"=paste(tmp_df[tmp_df$name == "group..attrs.z","value"]),
    				 "expect"=paste(tmp_df[tmp_df$name == "group..attrs.expect","value"]),
    				 "index_mgf"=paste(tmp_df[tmp_df$name == "group..attrs.id","value"]),
    				 "modification"=paste(paste(paste(tmp_df[tmp_df$name == "group.protein.peptide.domain.aa.type","value"]),paste(tmp_df[tmp_df$name == "group.protein.peptide.domain.aa.at","value"]),paste(tmp_df[tmp_df$name == "group.protein.peptide.domain.aa.modified","value"]),sep=","),collapse=";"),
    				 "start"=unique(paste(tmp_df[tmp_df$name == "group.protein.peptide.domain..attrs.start","value"])),
    				 "end"=unique(paste(tmp_df[tmp_df$name == "group.protein.peptide.domain..attrs.end","value"])),
    				 "hyperscore"=unique(paste(tmp_df[tmp_df$name == "group.protein.peptide.domain..attrs.hyperscore","value"])),
                     "missed_cleavages"=paste(tmp_df[tmp_df$name=="group.protein.peptide.domain..attrs.missed_cleavages","value"]))

                     df$scan<-str_replace_all(df$scan,"\n","")
                     df$scan<-str_replace_all(df$scan," ","_")


                     xml_data_df <- rbind(xml_data_df,df)


        }else if((length(paste(tmp_df[tmp_df$name == "group.protein.note.text","value"])) > 1) & (length(unique(paste(tmp_df[tmp_df$name == "group.protein.peptide.domain..attrs.seq","value"]))) == 1) &	(length(unique(paste(tmp_df[tmp_df$name == "group.protein.peptide.domain..attrs.start","value"]))) > 1) &
		(length(unique(paste(tmp_df[tmp_df$name == "group.protein.peptide.domain..attrs.end","value"]))) > 1)){

            for (l in 1:length(paste(tmp_df[tmp_df$name == "group.protein.note.text","value"]))){

                df<- data.frame("Protein"=paste(tmp_df[tmp_df$name == "group.protein.note.text","value"][l]),
					"PeptideSeq"=unique(paste(tmp_df[tmp_df$name == "group.protein.peptide.domain..attrs.seq","value"])),
					"scan"=paste(tmp_df[tmp_df$name == "group.group.note.text","value"]),
					"PrecursorMZ"=paste(tmp_df[tmp_df$name == "group..attrs.mh","value"]),
					"charge"=paste(tmp_df[tmp_df$name == "group..attrs.z","value"]),
					"expect"=paste(tmp_df[tmp_df$name == "group..attrs.expect","value"]),
					"index_mgf"=paste(tmp_df[tmp_df$name == "group..attrs.id","value"]),
					"modification"=paste(paste(paste(tmp_df[tmp_df$name == "group.protein.peptide.domain.aa.type","value"]),paste(tmp_df[tmp_df$name == "group.protein.peptide.domain.aa.at","value"]),paste(tmp_df[tmp_df$name == "group.protein.peptide.domain.aa.modified","value"]),sep=","),collapse=";"),
					"start"=paste(tmp_df[tmp_df$name == "group.protein.peptide.domain..attrs.start","value"][l]),
					"end"=paste(tmp_df[tmp_df$name == "group.protein.peptide.domain..attrs.end","value"][l]),
					"hyperscore"=unique(paste(tmp_df[tmp_df$name == "group.protein.peptide.domain..attrs.hyperscore","value"])),
                    "missed_cleavages"=paste(tmp_df[tmp_df$name=="group.protein.peptide.domain..attrs.missed_cleavages","value"][l]))

                    df$scan<-str_replace_all(df$scan,"\n","")
                    df$scan<-str_replace_all(df$scan," ","_")


                    xml_data_df <- rbind(xml_data_df,df)


            }


        }else if((length(paste(tmp_df[tmp_df$name == "group.protein.note.text","value"])) == 1) & (length(unique(paste(tmp_df[tmp_df$name == "group.protein.peptide.domain..attrs.seq","value"]))) > 1) &	(length(unique(paste(tmp_df[tmp_df$name == "group.protein.peptide.domain..attrs.start","value"]))) > 1) &
		(length(unique(paste(tmp_df[tmp_df$name == "group.protein.peptide.domain..attrs.end","value"]))) > 1)){
            #    si tengo un protein id pero más de una secuencia
            df<-data.frame("Protein"=as.character(),
            "PeptideSeq"=as.character(),
            "scan"=as.character(),
            "PrecursorMZ"=as.character(),
            "charge"=as.character(),
            "expect"=as.character(),
            "index_mgf"=as.character(),
            "modification"=as.character(),
            "start"=as.character(),
            "end"=as.character(),
            "hyperscore"=as.character(),
            "missed_cleavages"=as.character())

            for (j in 1:length(paste(tmp_df[tmp_df$name == "group.protein.peptide.domain..attrs.end","value"]))){

                df_tmp<- data.frame("Protein"=paste(tmp_df[tmp_df$name == "group.protein.note.text","value"]),
					"PeptideSeq"=paste(tmp_df[tmp_df$name == "group.protein.peptide.domain..attrs.seq","value"][j]),
					"scan"=paste(tmp_df[tmp_df$name == "group.group.note.text","value"]),
					"PrecursorMZ"=paste(tmp_df[tmp_df$name == "group..attrs.mh","value"]),
					"charge"=paste(tmp_df[tmp_df$name == "group..attrs.z","value"]),
					"expect"=paste(tmp_df[tmp_df$name == "group..attrs.expect","value"]),
					"index_mgf"=paste(tmp_df[tmp_df$name == "group..attrs.id","value"]),
					"modification"=paste(paste(paste(tmp_df[tmp_df$name == "group.protein.peptide.domain.aa.type","value"]),paste(tmp_df[tmp_df$name == "group.protein.peptide.domain.aa.at","value"]),paste(tmp_df[tmp_df$name == "group.protein.peptide.domain.aa.modified","value"]),sep=","),collapse=";"),
					"start"=paste(paste(tmp_df[tmp_df$name == "group.protein.peptide.domain..attrs.start","value"][j])),
					"end"=paste(paste(tmp_df[tmp_df$name == "group.protein.peptide.domain..attrs.end","value"][j])),
					"hyperscore"=unique(paste(tmp_df[tmp_df$name == "group.protein.peptide.domain..attrs.hyperscore","value"])),
                    "missed_cleavages"=paste(paste(tmp_df[tmp_df$name=="group.protein.peptide.domain..attrs.missed_cleavages","value"][j])))

                    df<-rbind(df,df_tmp)
                    df<-unique(df)

                    df$scan<-str_replace_all(df$scan,"\n","")
                    df$scan<-str_replace_all(df$scan," ","_")

            }


        }else if((length(paste(tmp_df[tmp_df$name == "group.protein.note.text","value"])) > 1) & (length(unique(paste(tmp_df[tmp_df$name == "group.protein.peptide.domain..attrs.seq","value"]))) == 1) & (length(unique(paste(tmp_df[tmp_df$name == "group.protein.peptide.domain..attrs.start","value"]))) == 1) & (length(unique(paste(tmp_df[tmp_df$name == "group.protein.peptide.domain..attrs.end"
		,"value"]))) == 1)){

            #si tengo mas de un protein id y una secuencia
            df<-data.frame()

            for (j in 1:length(paste(tmp_df[tmp_df$name == "group.protein.peptide.domain..attrs.seq","value"]))){

                df_tmp<- data.frame("Protein"=paste(paste(tmp_df[tmp_df$name == "group.protein.note.text","value"][j])),
					"PeptideSeq"=unique(paste(tmp_df[tmp_df$name == "group.protein.peptide.domain..attrs.seq","value"])),
					"scan"=paste(tmp_df[tmp_df$name == "group.group.note.text","value"]),
					"PrecursorMZ"=paste(tmp_df[tmp_df$name == "group..attrs.mh","value"]),
					"charge"=paste(tmp_df[tmp_df$name == "group..attrs.z","value"]),
					"expect"=paste(tmp_df[tmp_df$name == "group..attrs.expect","value"]),
					"index_mgf"=paste(tmp_df[tmp_df$name == "group..attrs.id","value"]),
					"modification"=paste(paste(paste(tmp_df[tmp_df$name == "group.protein.peptide.domain.aa.type","value"]),paste(tmp_df[tmp_df$name == "group.protein.peptide.domain.aa.at","value"]),paste(tmp_df[tmp_df$name == "group.protein.peptide.domain.aa.modified","value"]),sep=","),collapse=";"), "start"=unique(paste(tmp_df[tmp_df$name == "group.protein.peptide.domain..attrs.start","value"])),
					"end"=unique(paste(tmp_df[tmp_df$name == "group.protein.peptide.domain..attrs.end","value"])),
					"hyperscore"=unique(paste(tmp_df[tmp_df$name == "group.protein.peptide.domain..attrs.hyperscore","value"])),
                    "missed_cleavages"=paste(tmp_df[tmp_df$name=="group.protein.peptide.domain..attrs.missed_cleavages","value"][j]))

                    df<-rbind(df,df_tmp)
                    df<-unique(df)

                    df$scan<-str_replace_all(df$scan,"\n","")
                    df$scan<-str_replace_all(df$scan," ","_")

            }
            xml_data_df <- rbind(xml_data_df,df)


        }else if((length(paste(tmp_df[tmp_df$name == "group.protein.note.text","value"])) == 1) & (length(unique(paste(tmp_df[tmp_df$name == "group.protein.peptide.domain..attrs.seq","value"]))) == 1) & (length(unique(paste(tmp_df[tmp_df$name == "group.protein.peptide.domain..attrs.start","value"]))) > 1) & (length(unique(paste(tmp_df[tmp_df$name == "group.protein.peptide.domain..attrs.end"
		,"value"]))) > 1)){

            for(m in 1:length(paste(tmp_df[tmp_df$name == "group.protein.peptide.domain..attrs.start","value"]))){

                df<- data.frame("Protein"=paste(tmp_df[tmp_df$name == "group.protein.note.text","value"]),
					"PeptideSeq"=unique(paste(tmp_df[tmp_df$name == "group.protein.peptide.domain..attrs.seq","value"])),
					"scan"=paste(tmp_df[tmp_df$name == "group.group.note.text","value"]),
					"PrecursorMZ"=paste(tmp_df[tmp_df$name == "group..attrs.mh","value"]),
					"charge"=paste(tmp_df[tmp_df$name == "group..attrs.z","value"]),
					"expect"=paste(tmp_df[tmp_df$name == "group..attrs.expect","value"]),
					"index_mgf"=paste(tmp_df[tmp_df$name == "group..attrs.id","value"]),
					"modification"=paste(paste(paste(tmp_df[tmp_df$name == "group.protein.peptide.domain.aa.type","value"]),paste(tmp_df[tmp_df$name == "group.protein.peptide.domain.aa.at","value"]),paste(tmp_df[tmp_df$name == "group.protein.peptide.domain.aa.modified","value"]),sep=","),collapse=";"),
                    "start"=paste(unique(paste(tmp_df[tmp_df$name == "group.protein.peptide.domain..attrs.start","value"][m]))),
					"end"=paste(unique(paste(tmp_df[tmp_df$name == "group.protein.peptide.domain..attrs.end","value"][m]))),
					"hyperscore"=paste(unique(paste(tmp_df[tmp_df$name == "group.protein.peptide.domain..attrs.hyperscore","value"][m]))),
                    "missed_cleavages"=paste(tmp_df[tmp_df$name=="group.protein.peptide.domain..attrs.missed_cleavages","value"][m]))

                    df$scan<-str_replace_all(df$scan,"\n","")
                    df$scan<-str_replace_all(df$scan," ","_")


                xml_data_df <- rbind(xml_data_df,df)


            }

        }else if((length(paste(tmp_df[tmp_df$name == "group.protein.note.text","value"])) > 1) &	(length(unique(paste(tmp_df[tmp_df$name == "group.protein.peptide.domain..attrs.seq","value"]))) > 1)){
        # si tengo más de un protein id y más de una secuencia


            if((length(paste(tmp_df[tmp_df$name == "group.protein.note.text","value"]))) == (length(paste(tmp_df[tmp_df$name == "group.protein.peptide.domain..attrs.seq","value"])))){

                for (j in 1:length(paste(tmp_df[tmp_df$name == "group.protein.peptide.domain..attrs.seq","value"]))){

                    df<- data.frame("Protein"=paste(tmp_df[tmp_df$name == "group.protein.note.text","value"][j]),
						"PeptideSeq"=paste(tmp_df[tmp_df$name == "group.protein.peptide.domain..attrs.seq","value"][j]),
						"scan"=paste(tmp_df[tmp_df$name == "group.group.note.text","value"]),
						"PrecursorMZ"=paste(tmp_df[tmp_df$name == "group..attrs.mh","value"]),
						"charge"=paste(tmp_df[tmp_df$name == "group..attrs.z","value"]),
						"expect"=paste(tmp_df[tmp_df$name == "group..attrs.expect","value"]),
						"index_mgf"=paste(tmp_df[tmp_df$name == "group..attrs.id","value"]),
						"modification"=paste(paste(paste(tmp_df[tmp_df$name == "group.protein.peptide.domain.aa.type","value"]),paste(tmp_df[tmp_df$name == "group.protein.peptide.domain.aa.at","value"]),paste(tmp_df[tmp_df$name == "group.protein.peptide.domain.aa.modified","value"]),sep=","),collapse=";"),
                        "start"=paste(paste(tmp_df[tmp_df$name == "group.protein.peptide.domain..attrs.start","value"][j])),
						"end"=paste(paste(tmp_df[tmp_df$name == "group.protein.peptide.domain..attrs.end","value"][j])),
						"hyperscore"=paste(unique(paste(tmp_df[tmp_df$name == "group.protein.peptide.domain..attrs.hyperscore","value"][j]))),
                        "missed_cleavages"=paste(paste(tmp_df[tmp_df$name=="group.protein.peptide.domain..attrs.missed_cleavages","value"][j])))
                    #df$scan<-str_replace(df$scan,"\n","")
                    df$scan<-str_replace_all(df$scan,"\n","")
                    df$scan<-str_replace_all(df$scan," ","_")

                        xml_data_df <- rbind(xml_data_df,df)

                }



            }else{


                df<-data.frame()
                for (k in 1:length(paste(tmp_df[tmp_df$name == "group.protein.peptide.domain..attrs.seq","value"]))){

                    df_temp<- data.frame("Protein"=paste(tmp_df[tmp_df$name == "group.protein.note.text","value"][as.numeric(unlist(strsplit(paste(tmp_df[tmp_df$name == "group.protein.peptide.domain..attrs.id","value"][k]),"[.]"))[2])]),
						"PeptideSeq"=paste(tmp_df[tmp_df$name == "group.protein.peptide.domain..attrs.seq","value"][k]),
						"scan"=paste(tmp_df[tmp_df$name == "group.group.note.text","value"]),
						"PrecursorMZ"=paste(tmp_df[tmp_df$name == "group..attrs.mh","value"]),
						"charge"=paste(tmp_df[tmp_df$name == "group..attrs.z","value"]),
						"expect"=paste(tmp_df[tmp_df$name == "group..attrs.expect","value"]),
						"index_mgf"=paste(tmp_df[tmp_df$name == "group..attrs.id","value"]),
						"modification"=paste(paste(paste(tmp_df[tmp_df$name == "group.protein.peptide.domain.aa.type","value"]),paste(tmp_df[tmp_df$name == "group.protein.peptide.domain.aa.at","value"]),paste(tmp_df[tmp_df$name == "group.protein.peptide.domain.aa.modified","value"]),sep=","),collapse=";"), "start"=paste(tmp_df[tmp_df$name == "group.protein.peptide.domain..attrs.start","value"][k]),
						"end"=paste(tmp_df[tmp_df$name == "group.protein.peptide.domain..attrs.end","value"][k]),
						"hyperscore"=unique(paste(tmp_df[tmp_df$name == "group.protein.peptide.domain..attrs.hyperscore","value"][k])),
                        "missed_cleavages"=paste(tmp_df[tmp_df$name=="group.protein.peptide.domain..attrs.missed_cleavages","value"][k]))

                        df_temp<-unique(df_temp)
                        df<-rbind(df,df_temp)
                }
                xml_data_df <- rbind(xml_data_df,df)

            }



        }



    }

    return(xml_data_df)
}

xtandem_output<-parsedataXML(xml_data_unlist)

if(nrow(xtandem_output)==0){
    a<-cat("Revisa el bucle de parseo","\n")
    write.table(a, file = paste(gsub(".xml","",xmlFile), ".fail_txt", sep = ""), quote = FALSE, row.names = FALSE, col.names = TRUE, sep="\t")

}

#save(xml_data_df, file = paste(gsub(".xml","",xmlFile), ".rda", sep = ""))
write.table(xtandem_output, file = tsv_out, quote = FALSE, row.names = FALSE, col.names = TRUE, sep="\t")
cat(tsv_out,"\n")
cat(system("date"))
