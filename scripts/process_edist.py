#process_edist.py
# argv[1]: fasta file containing proteome
# argv[2]: name and path of hmmerscan results
# argv[3]: level of hmms scanned against (ie eukaryotes)
# argv[4]: name and path for output (can make new file or add to existing)
# Returns a space seperated text file with level, proteinID, orthogroup1(top hit), evalue1, orthogroup2 (2nd hit), evalue2 and proteomeID (from proteome file name). If no hits orthogroup1 is listed as proteinID.
# For analysis of evalue distributions.

import os.path
import sys
from Bio import SearchIO

def process_edist():
	#get the name of the proteome from the file name
	pid = sys.argv[1].replace(".fasta", "")
	pid = pid.split('/')
	pid = pid[len(pid)-1] 
	#read in results
	results = SearchIO.parse(sys.argv[2], "hmmer3-text")
	#initialize list to add entries to
	processed = []

	for protein in results:
		if len(protein) == 0: #if a query has no hits groupid=proteinid
			OG1id = protein.id
			e1 = "n/a"
			OG2id = "n/a"
			e2 = "n/a"
		elif len(protein) == 1: #if a query has 1 hit og2 and eval2 are n/a
			OG1 = protein[0]
			OG1id = OG1.id.replace(".meta_raw", "")
			e1 = str(OG1.evalue)
			OG2id = "n/a"
			e2 = "n/a"
		else:	#if a query has more hits the top two are recorded
			OG1 = protein[0]
			OG1id = OG1.id.replace(".meta_raw", "")
			e1 = str(OG1.evalue)
			OG2 = protein[1]
			OG2id = OG2.id.replace(".meta_raw", "")
			e2 = str(OG2.evalue)
		processed.append((sys.argv[3]+" ", protein.id+" ", OG1id+" ", e1+" ", OG2id+" ", e2+" ", pid))	
	
	#write entries to file
	if os.path.isfile(sys.argv[4]): #if an existing file was provided append entries to that file
		output = open(sys.argv[4], "a")
	else: #else make a new file and add a header
		output = open(sys.argv[4], "w")
		output.write("Level ProteinID Group1 evalue1 Group2 evalue2 ProteomeID \n")
	for i in processed:
		output.write("".join(str(s) for s in i) + "\n")
	output.close()	

#check for correct number of inputs
if len(sys.argv) !=5:
	print("Incorrect inputs, needs proteome.fasta, path/to/hmmer results, level, and output name.")
else:
	process_edist()

