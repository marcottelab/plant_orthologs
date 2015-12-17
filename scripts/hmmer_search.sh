#!/bin/bash
#hmmer_search.sh
# $1: fasta file containing proteome
# $2: path and name of directory to store results 
# $3: name for hmm database, if already created include path
# optional $4: directory where hmms are stored if haven't made hmm database yet

if [[ $# -eq 4 ]]; then
	cat $4/*.hmm > $2/$3
	hmmpress $2/$3
	hmmsearch --tblout $2/"search_results.tsv" $2/$3 $1
else
	hmmsearch --tblout $2/"search_results.tsv" $3 $1
fi	
