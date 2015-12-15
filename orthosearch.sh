#!/bin/bash
# takes a group of hmms and a proteome
# returns a proteome with the sequences that map to the same hmm concatonated
# $1: directory where hmms are stored
# $2: fasta file containing proteome 
# $3: name of directory to save output to
# $4: name of fasta file that will be created

mkdir $3 
source hmmer_search.sh $1 $2 $3 
python process_search.py $2 $3 $4

 
