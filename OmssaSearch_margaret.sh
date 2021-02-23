#!/bin/bash


CURRENTFOLDER=$1	# /mnt/beegfs/agarin/dato-activo/03_Analysis/agarin/26_Navarrabiomed_Missing_Ene18/
DATASET=$2			# /mnt/beegfs/agarin/dato-activo/03_Analysis/agarin/26_Navarrabiomed_Missing_Ene18/MGFFiles/349-028
DBINDEXROOT=$3	 	# /mnt/beegfs/agarin/dato-activo/03_Analysis/agarin/26_Navarrabiomed_Missing_Ene18/uniprot_sprot_2017_12
QUALITY=$4 			# HIGH; MEDIUM OR LOW.

echo "CURRENTFOLDER --> "$CURRENTFOLDER
echo "DATASET --> "$DATASET
echo "DBINDEXROOT --> "$DBINDEXROOT
echo "QUALITY --> "$QUALITY

datasetname=$(basename "${DATASET}")
omssafolder=$CURRENTFOLDER'Omssa_files'
mgfFolder=$CURRENTFOLDER'MGFFiles/'$datasetname
mkdir -p $omssafolder

echo "datasetname --> "$datasetname
echo "omssafolder --> "$omssafolder
echo "mgfFolder --> "$mgfFolder

mgfFiles=$(find "${mgfFolder}" -maxdepth 2 -iname "*.mgf")
mgfFiles=( $mgfFiles )
echo "mgfFiles --> "$mgfFiles
datasetOmssaTargetFolder=$omssafolder'/'$datasetname'/'
echo "datasetOmssaTargetFolder --> "$datasetOmssaTargetFolder
mkdir -p $datasetOmssaTargetFolder
datasetOmssaDecoyFolder=$omssafolder'/'$datasetname'-D/'
echo "datasetOmssaDecoyFolder --> "$datasetOmssaDecoyFolder
mkdir -p $datasetOmssaDecoyFolder

nMGFFiles=${#mgfFiles[@]}
echo "nMGFFiles --> "$nMGFFiles

OmssaFilesTarget=$(find "${datasetOmssaTargetFolder}" -iname "*.xml")
OmssaFilesTarget=( $OmssaFilesTarget )
nOmssaFilesTarget=${#OmssaFilesTarget[@]}
echo "datasetOmssaTargetFolder --> "$datasetOmssaTargetFolder
echo "OmssaFilesTarget --> "$OmssaFilesTarget
echo "nOmssaFilesTarget --> "$nOmssaFilesTarget

OmssaFilesDecoy=$(find "${datasetOmssaDecoyFolder}" -iname "*.xml")
OmssaFilesDecoy=( $OmssaFilesDecoy )
nOmssaFilesDecoy=${#OmssaFilesDecoy[@]}
echo "OmssaFilesDecoy --> "$OmssaFilesDecoy
echo "nOmssaFilesDecoy --> "$nOmssaFilesDecoy



if [ $nMGFFiles != $nOmssaFilesTarget ]; then
	echo 'Omssa target DB search...'
	for ((i=0; i < ${#mgfFiles[@]}; i ++))
	do
		mgfile=${mgfFiles[$i]}
		mgffileroot="${mgfile%.*}"
    #echo "mgffileroot --> "$mgffileroot
		DBINDEXFILE=$DBINDEXROOT'_index'
    #echo "DBINDEXFILE --> " $DBINDEXFILE
		# e trypsin; -i b,y ions; -mf carbamidomethyl C static mod; -mv oxidized M variable mod;
    if  [[ $QUALITY == "HIGH" ]]; then
      echo "high quality" # precursor tolerance 10 ppm; fragment tolerance 0.05; -hl retain top 10 hits; -v 1 cleavage
      omssacl -e 0 -i 1,4 -mf 3 -mv 1 -te 10 -teppm -to 0.05 -hl 10 -v 1 -fm $mgffileroot".mgf" -d $DBINDEXFILE -op $mgffileroot".pepXML" -ox $mgffileroot".xml" -oc $mgffileroot"_summary.csv"
    elif [[ $QUALITY == "MEDIUM" ]]; then
      echo "medium quality" # precursor tolerance 20 ppm; fragment tolerance 0.05; -hl retain top 10 hits; -v 1 cleavage
      omssacl -e 0 -i 1,4 -mf 3 -mv 1 -te 20 -teppm -to 0.05 -hl 10 -v 1 -fm $mgffileroot".mgf" -d $DBINDEXFILE -op $mgffileroot".pepXML" -ox $mgffileroot".xml" -oc $mgffileroot"_summary.csv"
    else
      echo "low quality" # precursor tolerance 20 ppm; fragment tolerance 0.5; -hl retain top 10 hits; -v 1 cleavage
      omssacl -e 0 -i 1,4 -mf 3 -mv 1 -te 20 -teppm -to 0.5 -hl 10 -v 1 -fm $mgffileroot".mgf" -d $DBINDEXFILE -op $mgffileroot".pepXML" -ox $mgffileroot".xml" -oc $mgffileroot"_summary.csv"
		fi
		mv "$mgffileroot.pepXML" $datasetOmssaTargetFolder
		mv "$mgffileroot.xml" $datasetOmssaTargetFolder
		mv $mgffileroot"_summary.csv" $datasetOmssaTargetFolder
	done
else
	echo "Target searches were already done. Continue"
fi

if [ $nMGFFiles != $nOmssaFilesDecoy ]; then
	echo 'Omssa decoy DB search...'
	for ((i=0; i < ${#mgfFiles[@]}; i ++))
	do
		mgfile=${mgfFiles[$i]}
		mgffileroot="${mgfile%.*}"
		#echo "mgfile --> "$mgfile
    #echo "mgffileroot --> "$mgffileroot
		DBINDEXFILE=$DBINDEXROOT'_D_index'
    #echo "DBINDEXFILE --> "$DBINDEXFILE
    # e trypsin; -i b,y ions; -mf carbamidomethyl C static mod; -mv oxidized M variable mod;
    if  [[ $QUALITY == "HIGH" ]]; then
      echo "high quality" # precursor tolerance 10 ppm; fragment tolerance 0.05; -hl retain top 10 hits; -v 1 cleavage
      omssacl -e 0 -i 1,4 -mf 3 -mv 1 -te 10 -teppm -to 0.05 -hl 10 -v 1 -fm $mgffileroot".mgf" -d $DBINDEXFILE -op $mgffileroot".pepXML" -ox $mgffileroot".xml" -oc $mgffileroot"_summary.csv"
    elif [[ $QUALITY == "MEDIUM" ]]; then
      echo "medium quality" # precursor tolerance 20 ppm; fragment tolerance 0.05; -hl retain top 10 hits; -v 1 cleavage
      omssacl -e 0 -i 1,4 -mf 3 -mv 1 -te 20 -teppm -to 0.05 -hl 10 -v 1 -fm $mgffileroot".mgf" -d $DBINDEXFILE -op $mgffileroot".pepXML" -ox $mgffileroot".xml" -oc $mgffileroot"_summary.csv"
    else
      echo "low quality" # precursor tolerance 20 ppm; fragment tolerance 0.5; -hl retain top 10 hits; -v 1 cleavage
      omssacl -e 0 -i 1,4 -mf 3 -mv 1 -te 20 -teppm -to 0.5 -hl 10 -v 1 -fm $mgffileroot".mgf" -d $DBINDEXFILE -op $mgffileroot".pepXML" -ox $mgffileroot".xml" -oc $mgffileroot"_summary.csv"
    fi
		#echo "mgffileroot.pepXML --> "$mgffileroot".pepXML"
		#echo "mgfileroot.xml --> "$mgffileroot".xml"
		#echo "mgffileroot_summary.csv --> "$mgffileroot"_summary.csv"
		#echo "datasetOmssaTargetFolder --> "$datasetOmssaDecoyFolder
		mv $mgffileroot".pepXML" $datasetOmssaDecoyFolder
		mv $mgffileroot".xml" $datasetOmssaDecoyFolder
		mv $mgffileroot"_summary.csv" $datasetOmssaDecoyFolder
	done
else
	echo "Decoy searches were already done. Continue"
fi
wait
