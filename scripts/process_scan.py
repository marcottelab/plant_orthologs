#process_scan.py 
#takes in hmmscan output table and returns text file with level, proteinID
#ortho group,and proteome id.
# argv[1]: fasta file containing proteome
# argv[2]: path to hmmer_scan.sh results
# argv[3]: level
# argv[4]: path and name for new text file or existing file to add to

import os.path
import sys
from Bio import SeqIO
from Bio import SearchIO

def process_scan():
	pid = sys.argv[1].replace(".fasta", "")
	pid = pid.split('/')
	pid = pid[len(pid)-1] 
	results = SearchIO.parse(sys.argv[2]+"/scan_results.txt", "hmmer3-text")
	processed = []

	for protein in results:
		if len(protein) == 0:
			OGid = protein.id
			e = "n/a"
			e2 = "n/a"
		elif len(protein) == 1: 
			OG = protein[0]
			OGid = OG.id.replace(".meta_raw", "")
			e = str(OG.evalue)
			e2 = "n/a"
		else:	
			OG = protein[0]
			OGid = OG.id.replace(".meta_raw", "")
			e = str(OG.evalue)
			OG2 = protein[1]
			e2 = str(OG2.evalue)
		processed.append((sys.argv[3]+" ", protein.id+" ", OGid+" ", e+" ", e2+" ", pid))	
	if os.path.isfile(sys.argv[4]):
		output = open(sys.argv[4], "a")
	else:
		output = open(sys.argv[4], "w")
		output.write("Level ProteinID GroupID evalue evalue2 ProteomeID \n")
		
	for i in processed:
		output.write("".join(str(s) for s in i) + "\n")
	output.close()	

if len(sys.argv) !=5:
	print("Incorrect inputs, needs proteome.fasta, path/to/hmmer results, level, and output name.")
else:
	process_scan()
