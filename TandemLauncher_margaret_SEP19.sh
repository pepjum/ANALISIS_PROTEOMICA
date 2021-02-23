#! /bin/bash

# ./TandemLaucher_.sh /home/nostromo/data//proteomica/ejemplo_experimento/ uniprot_sprot_2017_12_CRAP.fasta uniprot_sprot_2017_12_CRAP_D.fasta 3 DECOY /home/nostromo/data//proteomica/ejemplo_experimento/MGFFiles/ HIGH

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
		EXP=$MGFFOLDER'/mgfFolderList.txt'
		# echo $EXP
		# echo $MGFFOLDER > $EXP
		find $MGFFOLDER -mindepth 1 -maxdepth 1 -type d -print > $EXP
      #
		dbtarget=$currentfolder$dbtargetfile
		dbtargetname="${dbtargetfile%.*}"
		dbdecoy=$currentfolder$dbdecoyfile
		dbdecoyname="${dbdecoyfile%.*}"
      #
		tandemFolder=$currentfolder'Tandem_Files/'
		mkdir -p $tandemFolder
		# PREPARE THE GENERAL PARAMETER FILES
		defaulinputFile="/opt/tandem-linux-17-02-01-4/bin/default_input_"$QUALITY".xml"
		defaulinputFileD="/opt/tandem-linux-17-02-01-4/bin/default_input_DECOY_"$QUALITY".xml"
		cp $defaulinputFile $tandemFolder
		cp $defaulinputFileD $tandemFolder
		replace "databaseURL" $dbtarget < /opt/tandem-linux-17-02-01-4/bin/taxonomy.xml > $tandemFolder"taxonomy.xml"
		replace "decoydatabaseURL" $dbdecoy < /opt/tandem-linux-17-02-01-4/bin/taxonomy_DECOY.xml > $tandemFolder"taxonomy_DECOY.xml"
		/home/margaret/data/pepe/scripts/ANALISIS_PROTEOMICA/PeptideMatch.sh "${currentfolder}" "${dbtargetfile}" "${dbdecoyfile}" "${decoyid}"
      #
		echo "EMPIEZA SEARCH"
		while IFS='' read -r line || [[ -n "$line" ]]; do
		  dataset=$line
		  echo "dataset -> "$dataset
		  datasetname=$(basename "${dataset}")
		  echo "datasetname -> "$datasetname
		  dbrootpath=$currentfolder$dbtargetname
		  echo "dbrootpath -> "$dbrootpath
		  txtFileName=$tandemFolder"Shotgun_tandem_"$datasetname".txt"
		  echo "txtFileName -> "$txtFileName
		  inputTargetFileName=$tandemFolder"input_"$datasetname".xml"
		  inputDecoyFileName=$tandemFolder"input_DECOY_"$datasetname".xml"
			echo "inputTargetFileName -> "$inputTargetFileName
			echo "inputDecoyFileName -> "$inputDecoyFileName
		  replace "default_input.xml" $tandemFolder"default_input_"$QUALITY".xml" < /opt/tandem-linux-17-02-01-4/bin/input.xml > $inputTargetFileName".tmp1"
		  replace "taxonomy.xml" $tandemFolder"taxonomy.xml" < $inputTargetFileName".tmp1" > $inputTargetFileName".tmp"

		  replace "defaultdecoyinputURL" $tandemFolder"default_input_DECOY_"$QUALITY".xml" < /opt/tandem-linux-17-02-01-4/bin/input_DECOY.xml > 	$inputDecoyFileName".tmp1"
		  replace "decoytaxonomyURL" $tandemFolder"taxonomy_DECOY.xml" < $inputDecoyFileName".tmp1" > $inputDecoyFileName".tmp"
		  echo "${currentfolder}" "${dataset}" "${inputTargetFileName}" "${inputDecoyFileName}"
		  /home/margaret/data/pepe/scripts/ANALISIS_PROTEOMICA/TandemSearch_margaret.sh "${currentfolder}" "${dataset}" "${inputTargetFileName}" "${inputDecoyFileName}"
		done < "${EXP}"
		echo "TERMINA SEARCH"

		echo "EMPIEZA PARSE"
		while IFS='' read -r line || [[ -n "$line" ]]; do
			dataset=$line
			datasetname=$(basename "${dataset}")
			echo "${currentfolder}" "${datasetname}"
			/home/margaret/data/pepe/scripts/ANALISIS_PROTEOMICA/TandemParser_margaret.sh "${currentfolder}" "${datasetname}"
		done < "${EXP}"
		echo "TERMINA PARSE"
		echo "EMPIEZA PEP TO PROTEIN"
		while IFS='' read -r line || [[ -n "$line" ]]; do
			dataset=$line
			datasetname=$(basename "${dataset}")
			txtFileName=$tandemFolder"Shotgun_tandem_"$datasetname".txt"

			/home/margaret/data/pepe/scripts/ANALISIS_PROTEOMICA/TandemPep2Prot_margaret.sh "${currentfolder}" "${dataset}" "${txtFileName}" "${dbtargetfile}" "${decoyid}"
		done < "${EXP}"
	  echo "TERMINA PEP TO PROTEIN"
		echo "EMPIEZA SHOTGUN ANALYSIS"
		while IFS='' read -r line || [[ -n "$line" ]]; do
		  dataset=$line
		  datasetname=$(basename "${dataset}")
		  txtFileName=$tandemFolder"Shotgun_tandem_"$datasetname".txt"
			echo "Rscript /home/margaret/data/pepe/scripts/ANALISIS_PROTEOMICA/ShotgunAnalysisAutomated_margaret.R "${txtFileName}" "${datasetname}" 0 "${decoyid}" 3 3 9 30"
			Rscript /home/margaret/data/pepe/scripts/ANALISIS_PROTEOMICA/ShotgunAnalysisAutomated_margaret.R "${txtFileName}" "${datasetname}" 0 "${decoyid}" 3 3 9 50 &
		done < "${EXP}"
		echo "TERMINA SHOTGUN ANALYSIS"
fi
