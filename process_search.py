#process_search.py
#take in results of hmmsearch.sh and 
#return fasta with sequences that map to same hmm concatonated
# argv[1]: fasta file containing proteome 
# argv[2]: directory where output of hmmsearch.sh is stored
# argv[3]: name of fasta file that will be created

import sys
from Bio import SearchIO
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

results = SearchIO.parse(sys.argv[2]+"/search_results.tsv", "hmmer3-tab")
proteins = SeqIO.to_dict(SeqIO.parse(sys.argv[1], "fasta")) 
records = []

for result in results:
	catseq = ""
	catname = ""
    	for hit in result:       	
		catseq += proteins[hit.id].seq
		catname += hit.id
	records.append(SeqRecord(catseq, id = result.id, description = catname))

SeqIO.write(records, sys.argv[3], "fasta")  

