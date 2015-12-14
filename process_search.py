#process_search.py
#process results from hmmsearch

from Bio import SearchIO
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet import generic_protein

results = SearchIO.parse("search_results.tsv", "hmmer3-tab")
proteome = SeqIO.parse("miniproteome.fasta", "fasta") 
records = []

for result in results:
	catseq = ""
	catname = ""
	for hit in result:
		hitseq = ""
		for protein in proteome:
			if hit.id == protein.id:
				hitseq += protein.seq
			else:
				continue
		catseq += hitseq
		catname += hit.id
	records.append(SeqRecord(catseq,
		 id = result.id, description = catname))

SeqIO.write(records, "catproteome.fasta", "fasta")  

