#!/bin/bash

#./TandemParser_mlsanchez.sh /home/nostromo/data/mlsanchez/proteomica/ejemplo_experimento/ Exp1_M5
#### ./TandemParser_mlsanchez.sh /home/margaret/data/mlsanchez/ejemplo_experimento/ Exp1_M5


# cambiar ruta en module.sh
CURRENTFOLDER=$1	# /mnt/beegfs/agarin/dato-activo/03_Analysis/agarin/26_Navarrabiomed_Missing_Ene18/
DATASETNAME=$2		# 349-028

tandemTargetFolder=$CURRENTFOLDER'Tandem_Files/'$DATASETNAME'/'
xmlFilesT=$(find "${tandemTargetFolder}" -maxdepth 1 -iname "*.xml")
echo $xmlFilesT
xmlFilesT=( $xmlFilesT )

for ((i=0; i < ${#xmlFilesT[@]}; i ++))
do
		xmlFilesTfile=${xmlFilesT[$i]}
		xmlFilesTfilenoext="${xmlFilesTfile%.*}"
		#echo $xmlFilesTfile -o $xmlFilesTfilenoext".tsv"
		# /opt/Lego/module.sh Parser -i $xmlFilesTfile -o $xmlFilesTfilenoext".tsv"
		Rscript /home/margaret/data/pepe/scripts/ANALISIS_PROTEOMICA/parseXML_XTANDEM_SEP19.R $xmlFilesTfile $xmlFilesTfilenoext".tsv"
done
wait


tandemDecoyFolder=$CURRENTFOLDER'Tandem_Files/'$DATASETNAME'-D/'
xmlFilesD=$(find "${tandemDecoyFolder}" -maxdepth 1 -iname "*.xml")
xmlFilesD=( $xmlFilesD )

for ((i=0; i < ${#xmlFilesD[@]}; i ++))
do
		xmlFilesDfile=${xmlFilesD[$i]}
		xmlFilesDfilenoext="${xmlFilesDfile%.*}"
		# /opt/Lego/module.sh Parser -i  $xmlFilesDfile -o $xmlFilesDfilenoext".tsv"
		Rscript /home/margaret/data/pepe/scripts/ANALISIS_PROTEOMICA/parseXML_XTANDEM_SEP19.R $xmlFilesDfile $xmlFilesDfilenoext".tsv"

done
wait
