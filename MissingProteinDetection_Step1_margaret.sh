#!/bin/bash


CURRENTFOLDER=$1	# /mnt/beegfs/agarin/dato-activo/03_Analysis/agarin/28_PFortes_Shotgun_lncRNA_Feb18/
FASTAFILE=$2		# /mnt/beegfs/agarin/dato-activo/03_Analysis/agarin/23_neXtprot_20170801_Nov2017/nextProtDB20170801.fasta
NEXTPROTFOLDER=$3	# /mnt/beegfs/agarin/dato-activo/03_Analysis/agarin/23_neXtprot_20170801_Nov2017/

FASTAFILEFOLDER=$(dirname "${FASTAFILE}")
FASTAFILENAME=$(basename "${FASTAFILE}")
FASTAFILENAMENOEXT="${FASTAFILENAME%.*}"
FASTAFILENAMEEXT="${FASTAFILENAME##*.}"

echo "starting the digestion of the database"
perl /opt/proteogest/proteogest.pl -i $FASTAFILE -d -a -c trypsin -R M -W 15.99 -G 1 -L 1

fileroot=$FASTAFILEFOLDER'/'$FASTAFILENAMENOEXT

echo "parsing proteogest"
Rscript /home/margaret/data/pepe/scripts/ANALISIS_PROTEOMICA/DigestDatabase.R $fileroot

echo "starting to process the database"

Rscript /home/margaret/data/pepe/scripts/ANALISIS_PROTEOMICA/ProcessNextprot_margaret.R $NEXTPROTFOLDER

echo "starting the detection of missing proteins"

Rscript /home/margaret/data/pepe/scripts/ANALISIS_PROTEOMICA/MissingProteinDetection_Step1_margaret.R $CURRENTFOLDER $fileroot $NEXTPROTFOLDER

echo "The detection of missing proteins has finished."
