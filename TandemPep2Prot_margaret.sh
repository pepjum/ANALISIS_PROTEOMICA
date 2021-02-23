#!/bin/bash

# ./TandemPep2Prot_mlsanchez.sh /home/nostromo/data/mlsanchez/proteomica/ejemplo_experimento /home/nostromo/data/mlsanchez/proteomica/ejemplo_experimento/MGFFiles/Exp1_M3 /home/nostromo/data/mlsanchez/proteomica/ejemplo_experimento/Tandem_Files/Shotgun_tandem_Exp1_M3.txt uniprot_sprot_2017_12_CRAP.fasta DECOY
# ./TandemPep2Prot_mlsanchez.sh /home/nostromo/data/mlsanchez/proteomica/EMBRIO_15_03_PXD003560-PXD006332/ /home/nostromo/data/mlsanchez/proteomica/EMBRIO_15_03_PXD003560-PXD006332/MGFFiles/PXD006314 /home/nostromo/data/mlsanchez/proteomica/EMBRIO_15_03_PXD003560-PXD006332/Tandem_Files/Shotgun_tandem_PXD006314.txt uniprot_sprot_2017_12_CRAP.fasta DECOY
# ./TandemPep2Prot_mlsanchez.sh /home/margaret/data/mlsanchez/28_PFortes_Shotgun_lncRNA_Feb18 /home/margaret/data/mlsanchez/28_PFortes_Shotgun_lncRNA_Feb18/MGFFiles/Exp1_M3 /home/margaret/data/mlsanchez/28_PFortes_Shotgun_lncRNA_Feb18/Tandem_Files/Shotgun_tandem_Exp1_M3.txt uniprot_sprot_2017_12_CRAP.fasta DECOY

CURRENTFOLDER=$1	# /mnt/beegfs/agarin/dato-activo/03_Analysis/agarin/26_Navarrabiomed_Missing_Ene18/
DATASET=$2			# /mnt/beegfs/agarin/dato-activo/03_Analysis/agarin/26_Navarrabiomed_Missing_Ene18/MGFFiles/349-028
FILENAME=$3 		# /mnt/beegfs/agarin/dato-activo/03_Analysis/agarin/26_Navarrabiomed_Missing_Ene18/Tandem_Files/Shotgun_tandem_349-028.txt
DBTARGETFILE=$4 	# uniprot_sprot_2017_12.fasta
DECOYID=$5			# DECOY

echo "CURRENTFOLDER --> "$CURRENTFOLDER
echo "DATASET --> "$DATASET
echo "FILENAME --> "$FILENAME
echo "DBTARGETFILE --> "$DBTARGETFILE
echo "DECOYID --> "$DECOYID

datasetname=$(basename "${DATASET}")
mgfFiles=$(find "${DATASET}" -maxdepth 2 -iname "*.mgf")
echo "mgfFiles --> "$mgfFiles
mgfFiles=( $mgfFiles )
nMGFFiles=${#mgfFiles[@]}

echo "datasetname --> "$datasetname
echo "nMGFFiles --> "$nMGFFiles

datasetTandemTargetFolder=$( echo $DATASET | sed -e  "s/MGFFiles/Tandem_Files/g")
datasetTandemDecoyFolder=$datasetTandemTargetFolder'-D/'
datasetTandemTargetFolder=$datasetTandemTargetFolder

echo "datasetTandemDecoyFolder --> "$datasetTandemDecoyFolder
echo "datasetTandemTargetFolder --> "$datasetTandemTargetFolder

TandemTargetFiles=$(find "${datasetTandemTargetFolder}" -iname "*.tsv"  -not -name "*_peptides.tsv" -not -name "*_peptides_out.tsv" -not -name "*_corrected.tsv" -not -name "*_peptides_log.log")
echo "TandemTargetFiles --> "$TandemTargetFiles
TandemTargetFiles=( $TandemTargetFiles )
nTandemTargetFiles=${#TandemTargetFiles[@]}
echo "nTandemTargetFiles --> "$nTandemTargetFiles

CorrectedTargetFiles=$(find "${datasetTandemTargetFolder}" -iname "*_corrected.tsv")
echo "CorrectedTargetFiles --> "$CorrectedTargetFiles
CorrectedTargetFiles=( $CorrectedTargetFiles )
nCorrectedTargetFiles=${#CorrectedTargetFiles[@]}
echo "nCorrectedTargetFiles --> "$nCorrectedTargetFiles


PeptideTargetFiles=$(find "${datasetTandemTargetFolder}" -iname "*_peptides.tsv")
echo "PeptideTargetFiles --> "$PeptideTargetFiles
PeptideTargetFiles=( $PeptideTargetFiles )
nPeptideTargetFiles=${#PeptideTargetFiles[@]}
echo "nPeptideTargetFiles --> "$nPeptideTargetFiles

PeptideOutTargetFiles=$(find "${datasetTandemTargetFolder}" -iname "*_peptides_out.tsv")
echo "PeptideOutTargetFiles --> "$PeptideOutTargetFiles
PeptideOutTargetFiles=( $PeptideOutTargetFiles )
nPeptideOutTargetFiles=${#PeptideOutTargetFiles[@]}
echo "nPeptideOutTargetFiles --> "$nPeptideOutTargetFiles

PeptideLogTargetFiles=$(find "${datasetTandemTargetFolder}" -iname "*_peptides_log.log")
echo "PeptideLogTargetFiles --> "$PeptideLogTargetFiles
PeptideLogTargetFiles=( $PeptideLogTargetFiles )
nPeptideLogTargetFiles=${#PeptideLogTargetFiles[@]}
echo "nPeptideLogTargetFiles --> "$nPeptideLogTargetFiles



TandemDecoyFiles=$(find "${datasetTandemDecoyFolder}" -iname "*.tsv"  -not -name "*_peptides.tsv" -not -name "*_peptides_out.tsv" -not -name "*_corrected.tsv" -not -name "*_peptides_log.log")
echo "TandemDecoyFiles --> "$TandemDecoyFiles
TandemDecoyFiles=( $TandemDecoyFiles )
nTandemDecoyFiles=${#TandemDecoyFiles[@]}
echo "nTandemDecoyFiles --> "$nTandemDecoyFiles

CorrectedDecoyFiles=$(find "${datasetTandemDecoyFolder}" -iname "*_corrected.tsv")
echo "CorrectedDecoyFiles --> "$CorrectedDecoyFiles
CorrectedDecoyFiles=( $CorrectedDecoyFiles )
nCorrectedDecoyFiles=${#CorrectedDecoyFiles[@]}
echo "nCorrectedDecoyFiles --> "$nCorrectedDecoyFiles

PeptideDecoyFiles=$(find "${datasetTandemDecoyFolder}" -iname "*_peptides.tsv")
echo "PeptideDecoyFiles --> "$PeptideDecoyFiles
PeptideDecoyFiles=( $PeptideDecoyFiles )
nPeptideDecoyFiles=${#PeptideDecoyFiles[@]}
echo "nPeptideDecoyFiles --> "$nPeptideDecoyFiles

PeptideOutDecoyFiles=$(find "${datasetTandemDecoyFolder}" -iname "*_peptides_out.tsv")
echo "PeptideOutDecoyFiles --> "$PeptideOutDecoyFiles
PeptideOutDecoyFiles=( $PeptideOutDecoyFiles )
nPeptideOutDecoyFiles=${#PeptideOutDecoyFiles[@]}
echo "nPeptideOutDecoyFiles --> "$nPeptideOutDecoyFiles

PeptideLogDecoyFiles=$(find "${datasetTandemDecoyFolder}" -iname "*_peptides_log.log")
echo "PeptideLogDecoyFiles --> "$PeptideLogDecoyFiles
PeptideLogDecoyFiles=( $PeptideLogDecoyFiles )
nPeptideLogDecoyFiles=${#PeptideLogDecoyFiles[@]}
echo "nPeptideLogDecoyFiles --> "$nPeptideLogDecoyFiles

target=0
if [ $nMGFFiles != $nTandemTargetFiles ]; then
	echo "There are no search files in the folder. Try agarin"
	exit
elif [ $nMGFFiles != $nCorrectedTargetFiles ]; then
	target=1
elif [ $nMGFFiles != $nPeptideTargetFiles ]; then
	target=1
elif [ $nMGFFiles != $nPeptideOutTargetFiles ]; then
	target=1
elif [ $nMGFFiles != $nPeptideLogTargetFiles ]; then
	target=1
else
	target=0
fi

decoy=0
if [ $nMGFFiles != $nTandemDecoyFiles ]; then
	echo "There are no search files in the folder. Try agarin"
	exit
elif [ $nMGFFiles != $nCorrectedDecoyFiles ]; then
	decoy=1
elif [ $nMGFFiles != $nPeptideDecoyFiles ]; then
	decoy=1
elif [ $nMGFFiles != $nPeptideOutDecoyFiles ]; then
	decoy=1
elif [ $nMGFFiles != $nPeptideLogDecoyFiles ]; then
	decoy=1
else
	decoy=0
fi

if [ $target -eq 1 ]; then
	if [ $decoy -eq 1 ]; then
		TOD=0
		echo "doing PeptideMatch for target and decoy files."
	else
		TOD=1
		echo "doing PeptideMatch for target files. Decoy files are already processed."
	fi
else
	if [ $decoy -eq 1 ]; then
		TOD=2
		echo "doing PeptideMatch for decoy files. Target files are already processed."
	else
		TOD=3
	fi
fi

if [ $TOD -ne 3 ]; then
	# TOD; 1 si solo hay que hacer target; 2 si solo hay que hacer decoy; 0 si hay que hacer los dos.
	Rscript /home/margaret/data/pepe/scripts/ANALISIS_PROTEOMICA/TandemPep2Prot_margaret.R $CURRENTFOLDER $datasetname $FILENAME $DBTARGETFILE $DECOYID $TOD
else
	echo "PeptideMatch is already done. Continue."
fi
