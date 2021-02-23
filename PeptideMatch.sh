#!/bin/bash

# ./PeptideMatch_mlsanchez /home/nostromo/data/mlsanchez/proteomica/ejemplo_experimento/ uniprot_sprot_2017_12_CRAP.fasta uniprot_sprot_2017_12_CRAP_D.fasta DECOY

### ./PeptideMatch_mlsanchez /home/margaret/data/mlsanchez/ejemplo_experimento/ uniprot_sprot_2017_12_CRAP.fasta uniprot_sprot_2017_12_CRAP_D.fasta DECOY


CURRENTFOLDER=$1 	# /mnt/beegfs/agarin/dato-activo/03_Analysis/agarin/26_Navarrabiomed_Missing_Ene18/
DBTARGETFILE=$2 	# uniprot_sprot_2017_12.fasta
DBDECOYFILE=$3 		# uniprot_sprot_2017_12_D.fasta
DECOYID=$4			# DECOY

echo 'Currentfolder: '$CURRENTFOLDER
echo 'DBTargetFile: '$DBTARGETFILE
echo 'DBDecoyFile: '$DBDECOYFILE
echo 'DecoyID: '$DECOYID

dbname="${DBTARGETFILE%.*}"
peptidematchdir=$CURRENTFOLDER'PeptideMatch_'$dbname'_index'
dbtargetpath=$CURRENTFOLDER''$DBTARGETFILE

if [ ! -d "$peptidematchdir" ]; then
	echo 'doing target peptidematch'
	java -jar /opt/PeptideMatchCMD/PeptideMatchCMD_1.0.jar -a index -d $dbtargetpath -i $peptidematchdir -f
else
	echo 'target peptidematch already done. continue.'
fi

peptidematchdirD=$CURRENTFOLDER'PeptideMatch_'$dbname'_'$DECOYID'_index'
dbcometpath=$CURRENTFOLDER''$DBDECOYFILE

if [ ! -d "$peptidematchdirD" ]; then
	echo 'doing decoy peptidematch'
	java -jar /opt/PeptideMatchCMD/PeptideMatchCMD_1.0.jar -a index -d $dbcometpath -i $peptidematchdirD -f
else
	echo 'decoy peptidematch already done. continue.'
fi
wait
