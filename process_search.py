#process_search.py
#take in results of hmmsearch and return fasta with sequences of same OG contatonated

from Bio import SearchIO
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

results = SearchIO.parse("search_results.tsv", "hmmer3-tab")
proteins = SeqIO.to_dict(SeqIO.parse("miniproteome.fasta", "fasta")) 
records = []

for result in results:
	catseq = ""
	catname = ""
    	for hit in result:       	
		catseq += proteins[hit.id].seq
		catname += hit.id
	records.append(SeqRecord(catseq, id = result.id, description = catname))

SeqIO.write(records, "catproteome.fasta", "fasta")  

