#!/bin/bash

cat minihmm/*.hmm > minidata
hmmpress minidata
hmmscan --tblout  scan_results.tsv  minidata miniproteome.fasta

