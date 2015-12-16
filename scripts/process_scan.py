#process_scan.py 
#takes in hmmscan output table and returns .
from Bio import SeqIO
from Bio import SearchIO

results = SearchIO.parse("scan_results.tsv", "hmmer3-tab")

for result in results:
	

