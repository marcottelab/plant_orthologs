#!/bin/shell
#ortho_tabscan.sh 
# $1: fasta file containing proteome
# $2: path and name of directory to store results
# $3: name for hmm database, if already created include path
# $4: level
# $5: name and path for new text file
# optional $6: directory where hmms are stored if haven't made hmm database yet

if [[ $# -lt 5 ]] || [[ $3 -gt 6 ]]; then
	echo "You must enter 5 or 6 arguments"
else
	mkdir $2
	if [[ $# -eq 6 ]]; then
		source scripts/hmmer_scan.sh $1 $2 $3 $6
	else
		source scripts/hmmer_scan.sh $1 $2 $3
	fi
	python scripts/process_scan.py $1 $2 $4 $5
fi 
