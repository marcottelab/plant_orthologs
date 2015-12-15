#!/bin/bash
# $1: directory where hmms are stored
# $2: fasta file containing proteome 
# $3: directory to save output to

cat $1/*.hmm > $3/hmms
hmmpress $3/hmms
hmmsearch --tblout $3/search_results.tsv $3/hmms $2
