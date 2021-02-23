#!/bin/bash

#./TandemSearch_pepe.sh /home/nostromo/data/pepe/proteomica/ejemplo_experimento/ /home/nostromo/data/pepe/proteomica/ejemplo_experimento/MGFFiles/ /home/nostromo/data/pepe/proteomica/ejemplo_experimento/Tandem_Files/input_Exp1_M3.xml /home/nostromo/data/pepe/proteomica/ejemplo_experimento/Tandem_Files/input_DECOY_Exp1_M3.xml
#./TandemSearch_pepe.sh /home/nostromo/data/pepe/proteomica/EMBRIO_15_03_PXD003560-PXD006332/ /home/nostromo/data/pepe/proteomica/EMBRIO_15_03_PXD003560-PXD006332/MGFFiles/PXD001381 /home/nostromo/data/pepe/proteomica/EMBRIO_15_03_PXD003560-PXD006332/Tandem_Files/input_PXD001381.xml /home/nostromo/data/pepe/proteomica/EMBRIO_15_03_PXD003560-PXD006332/Tandem_Files/input_DECOY_PXD001381.xml
#./TandemSearch_pepe.sh /home/margaret/data/pepe/nci60/ /home/margaret/data/pepe/nci60/MGFFiles/00523_A10_P003818_B0A_A00_R1.mgf /home/margaret/data/pepe/nci60/Tandem_Files/input_MGFFiles.xml /home/margaret/data/pepe/nci60/Tandem_Files/input_DECOY_MGFFiles.xml



CURRENTFOLDER=$1			# /home/nostromo/data/pepe/script_proteomica/
DATASET=$2					# /home/nostromo/data/pepe/script_proteomica/MGFFiles/
INPUTTARGETFILENAME=$3	 	# /home/nostromo/data/pepe/script_proteomica/Tandem_Files/default_input_HIGH.xml
INPUTDECOYFILENAME=$4 		# /home/nostromo/data/pepe/script_proteomica/Tandem_Files/default_input_DECOY_HIGH.xml

echo "CURRENTFOLDER -> "$CURRENTFOLDER
echo "DATASET -> "$DATASET
echo "INPUTTARGETFILENAME -> "$INPUTTARGETFILENAME
echo "INPUTDECOYFILENAME -> "$INPUTDECOYFILENAME

datasetname=$(basename "${DATASET}")
echo "datasetname -> "$datasetname
tandemfolder=$CURRENTFOLDER'Tandem_Files/'
echo "tandemfolder -> "$tandemfolder
mgfFolder=$CURRENTFOLDER'MGFfiles/'$datasetname
echo "mgfFolder -> "$mgfFolder
#find $currentfolder$MGFFOLDER -mindepth 1 -maxdepth 1 -type f
mgfFiles=$(find "${DATASET}" -maxdepth 1 -type f -iname "*.mgf") # cambio 2 ya que los .mgf están en "un directorio más adentro"
echo "mgfFiles -> "$mgfFiles
mgfFiles=( $mgfFiles )


datasetTandemTargetFolder=$tandemfolder$datasetname
mkdir -p $datasetTandemTargetFolder
echo "datasetTandemTargetFolder -> "$datasetTandemTargetFolder
datasetTandemDecoyFolder=$tandemfolder$datasetname'-D/'
echo "datasetTandemDecoyFolder -> "$datasetTandemDecoyFolder
mkdir -p $datasetTandemDecoyFolder

nMGFFiles=${#mgfFiles[@]} #numeros de ficheros mgf que hay

TandemFilesTarget=$(find "${datasetTandemTargetFolder}" -iname "*.xml")
echo "TandemFilesTarget -> "$TandemFilesTarget
TandemFilesTarget=( $TandemFilesTarget )
nTandemFilesTarget=${#TandemFilesTarget[@]}
echo "nTandemFilesTarget -> "$nTandemFilesTarget
TandemFilesDecoy=$(find "${datasetTandemDecoyFolder}" -iname "*.xml")
echo "TandemFilesDecoy -> "$TandemFilesDecoy
TandemFilesDecoy=( $TandemFilesDecoy )
nTandemFilesDecoy=${#TandemFilesDecoy[@]}
echo "nTandemFilesDecoy -> "$nTandemFilesDecoy

if [ $nMGFFiles != $nTandemFilesTarget ]; then
	echo 'Tandem target DB search...'
	for ((i=0; i < ${#mgfFiles[@]}; i ++))
	do
		mgfile=${mgfFiles[$i]}
		mgffileroot="${mgfile%.*}"
		replace "test_spectra.mgf" $mgfile < $INPUTTARGETFILENAME".tmp" > $INPUTTARGETFILENAME".tmp1"
		replace "output.xml" $mgffileroot'.xml' < $INPUTTARGETFILENAME".tmp1" > $INPUTTARGETFILENAME
		/opt/tandem-linux-17-02-01-4/bin/tandem.exe $INPUTTARGETFILENAME
		xmlfilePath=$(dirname "${mgffileroot}")
		xmlFiles=$(find "${xmlfilePath}" -iname "*.xml")
		mv -f $xmlFiles $datasetTandemTargetFolder
	done

else
	echo "Target searches were already done. Continue"
fi

if [ -f $INPUTTARGETFILENAME".tmp" ]; then
   rm $INPUTTARGETFILENAME".tmp"
fi

if [ -f $INPUTTARGETFILENAME".tmp1" ]; then
   rm $INPUTTARGETFILENAME".tmp1"
fi

wait

if [ $nMGFFiles != $nTandemFilesDecoy ]; then
	echo 'Tandem decoy DB search...'
	for ((i=0; i < ${#mgfFiles[@]}; i ++))
	do
		mgfile=${mgfFiles[$i]}
		mgffileroot="${mgfile%.*}"
		replace "inputmgf" $mgfile < $INPUTDECOYFILENAME".tmp" > $INPUTDECOYFILENAME".tmp1"
		replace "decoyoutputxml" $mgffileroot'.xml' < $INPUTDECOYFILENAME".tmp1" > $INPUTDECOYFILENAME
		/opt/tandem-linux-17-02-01-4/bin/tandem.exe $INPUTDECOYFILENAME
		xmlfilePath=$(dirname "${mgffileroot}")
		xmlFiles_d=$(find  "${xmlfilePath}" -iname "*.xml")
		mv -f $xmlFiles_d $datasetTandemDecoyFolder
	done

else
	echo "Decoy searches were already done. Continue"
fi

if [ -f $INPUTDECOYFILENAME".tmp" ]; then
   rm -f $INPUTDECOYFILENAME".tmp"
fi

if [ -f $INPUTDECOYFILENAME".tmp1" ]; then
   rm -f $INPUTDECOYFILENAME".tmp1"
fi

wait
