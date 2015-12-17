#!/bin/bash
# ortho_catsearch.sh takes a group of hmms and a proteome
# returns a proteome with the sequences that map to the same hmm concatonated
# $1: fasta file containing proteome 
# $2: path and name for new directory to store results
# $3: name for hmm database, if already created include path
# $4: name and path for new fasta file
# optional $5: directory where hmms are stored if haven't made hmm database yet

if [[ $# -lt 4 ]] || [[ $# -gt 5 ]]; then
	echo 'You must enter 4 or 5 arguments'
else			
	mkdir $2
	if  [[ $# -eq 4 ]]; then 
		source scripts/hmmer_search.sh $1 $2 $3
	else   
		source scripts/hmmer_search.sh $1 $2 $3 $5
	fi
	python scripts/cat_search.py $1 $2 $4
fi 
