#proteome_breaker.py
#takes in fasta file returns a group of fasta files with a specified number of
#sequences each
# argv[1]: fasta file containing proteome
# argv[2]: number of sequences per output file
# argv[3]: directory to store output files

import sys
from Bio import SeqIO

def proteome_breaker():
	#read in proteome
	proteome = SeqIO.parse(sys.argv[1], "fasta")
	#get the name of the proteome
	pid = sys.argv[1].replace(".fasta", "")
	pid = pid.split('/')
	pid = pid[len(pid)-1] 
	
	n = 0
	entry = True

	while entry:
		batch = []	
		while len(batch) < int(sys.argv[2]):
			try:
				entry = proteome.next()
			except StopIteration:
				entry = None
			if entry is None:
				#End of file
				break
			batch.append(entry)
		if batch:
			SeqIO.write(batch, sys.argv[3]+"/"+pid+str(n)+".fasta", "fasta")
		n+=1	
	
#check for correct number of inputs
if len(sys.argv) !=4:
	print("Incorrect inputs, needs proteome.fasta, desired number of sequences per file, and directory to store output files.")
else:
	proteome_breaker()
