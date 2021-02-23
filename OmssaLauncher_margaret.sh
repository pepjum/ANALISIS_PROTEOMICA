#! /bin/bash

currentfolder=$1	# folder where all the data is stored. There must be one folder: MGFFiles
dbtargetfile=$2		# database file for target searches
dbdecoyfile=$3		# database file for decoy searches
dbrdafile=$4		# several options (rda file, 0, 1, 2, 3) => rda file: database rda file for MAYU protFDR analyses, 0: psmFDR, 1: protFDRMAYU, 2: pepFDR, 3: protFDR
decoyid=$5			# decoy sequence identifier.
MGFFOLDER=$6 		# mgfFolder with all the experiments to be done
QUALITY=$7			# HIGH; MEDIUM OR LOW.

echo "currentfolder --> "$currentfolder
echo "dbtargetfile --> "$dbtargetfile
echo "dbdecoyfile --> "$dbdecoyfile
echo "dbrdafile --> "$dbrdafile
echo "decoyid --> "$decoyid
echo "MGFFOLDER --> "$MGFFOLDER
echo "QUALITY --> "$QUALITY


if [[ $# == "7" ]]
	then
		EXP=$MGFFOLDER'mgfFolderList.txt'
		find $MGFFOLDER -mindepth 1 -maxdepth 1 -type d -print > $EXP

 		dbtarget=$currentfolder$dbtargetfile
		dbtargetname="${dbtargetfile%.*}"
 		dbdecoy=$currentfolder$dbdecoyfile
 		dbdecoyname="${dbdecoyfile%.*}"

		echo "dbtarget --> "$dbtarget
		echo "dbtargetname --> "$dbtargetname
		echo "dbdecoy --> "$dbdecoy
		echo "dbdecoyname --> "$dbdecoyname

		#PREPARE THE DATABASE
		blastdbtargetout=$currentfolder$dbtargetname'_index'
		/opt/ncbi-blast-2.2.29+/bin/makeblastdb -in $dbtarget -dbtype prot
		/opt/ncbi-blast-2.2.29+/bin/blastdb_aliastool -dblist $dbtarget -dbtype prot -title $dbtargetname -out $blastdbtargetout

		blastdbdecoyout=$currentfolder$dbdecoyname'_index'
		/opt/ncbi-blast-2.2.29+/bin/makeblastdb -in $dbdecoy -dbtype prot
		/opt/ncbi-blast-2.2.29+/bin/blastdb_aliastool -dblist $dbdecoy -dbtype prot -title $dbdecoyname -out $blastdbdecoyout

		echo "EMPIEZA SEARCH"

		while IFS='' read -r line || [[ -n "$line" ]]; do        # go through the file to get every line
		  dataset=$line
		  dbrootpath=$currentfolder$dbtargetname
			datasetname=$(basename "${dataset}")
		  # SEARCH
			echo "${currentfolder}" "${dataset}" "${dbrootpath}" "${QUALITY}"
		  /home/margaret/data/pepe/scripts/ANALISIS_PROTEOMICA/OmssaSearch_margaret.sh "${currentfolder}" "${datasetname}" "${dbrootpath}" "${QUALITY}"
		done < "${EXP}"

		echo "ACABA SEARCH"

		echo "EMPIEZA PARSE FROM XML TO TXT"

		while IFS='' read -r line || [[ -n "$line" ]]; do        # go through the file to get every line
		  dataset=$line
		  datasetname=$(basename "${dataset}")
		  txtFileName=$currentfolder"Omssa_files/Shotgun_omssa_"$datasetname".txt"
		  # PARSE FROM XML TO TXT
			echo "${currentfolder}" "${datasetname}"
		  /home/margaret/data/pepe/scripts/ANALISIS_PROTEOMICA/OmssaParser_margaret.sh "${currentfolder}" "${datasetname}"
		done < "${EXP}"

		echo "ACABA PARSE FROM XML TO TXT"

		echo "EMPIEZA CREATE THE TXT FILE"

		while IFS='' read -r line || [[ -n "$line" ]]; do        # go through the file to get every line
		  dataset=$line
		  datasetname=$(basename "${dataset}")
		  # CREATE THE TXT FILE
		  Rscript /home/margaret/data/pepe/scripts/ANALISIS_PROTEOMICA/OmssaCreateTxtFile_margaret.R "${currentfolder}" "${datasetname}"
		done < "${EXP}"

		#echo "ACABA CREATE THE TXT FILE"
		echo "EMPIEZA PEP TO PROTEIN"
		while IFS='' read -r line || [[ -n "$line" ]]; do
			dataset=$line
			datasetname=$(basename "${dataset}")
			txtFileName=$currentfolder"Omssa_files/Shotgun_omssa_"$datasetname".txt"

			/home/margaret/data/pepe/scripts/ANALISIS_PROTEOMICA/OmssaPep2Prot_margaret.sh "${currentfolder}" "${dataset}" "${txtFileName}" "${dbtargetfile}" "${decoyid}"
		done < "${EXP}"

		echo "EMPIEZA SHOTGUN ANALYSYS"

		while IFS='' read -r line || [[ -n "$line" ]]; do        # go through the file to get every line
		  dataset=$line
		  datasetname=$(basename "${dataset}")
		  txtFileName=$currentfolder"Omssa_files/Shotgun_omssa_"$datasetname".txt"
		  #echo ${txtFileName} ${datasetname}" 0 "${decoyid}" 3 4 9 30"
		  # SHOTGUN ANALYSIS
		  Rscript /home/margaret/data/pepe/scripts/ANALISIS_PROTEOMICA/ShotgunAnalysisAutomated_margaret.R "${txtFileName}" "${datasetname}" 0 "${decoyid}" 3 4 9 50

		done < "${EXP}"

		echo "ACABA SHOTGUN ANALYSYS"
fi
