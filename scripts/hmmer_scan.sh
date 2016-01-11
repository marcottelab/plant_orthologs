#!/bin/bash
#hmmer_scan.sh
# $1: fasta file containing proteome
# $2: path and name of directory to store results
# $3: name for hmm database, if already created include path
# $4: name for scan results file
# optional $5: directory where hmms are stored if haven't made hmm database yet

if [[ $# -eq 5 ]]; then
	echo $5/*.hmm | xargs cat > $2/$3
	hmmpress $2/$3
	hmmscan -o $2/$4 $2/$3 $1
else
	hmmscan -o $2/$4 $3 $1
fi
