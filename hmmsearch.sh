#!/bin/bash

cat minihmm/*.hmm > minidata
hmmpress minidata
hmmsearch --tblout search_results.tsv minidata miniproteome.fasta
