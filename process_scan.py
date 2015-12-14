#OGcat.py 
#takes in hmmscan output table and returns a fasta with sequences from the same OG concatonated.
from Bio import SeqIO
from Bio import SearchIO

proteome = list( SeqIO.parse("miniproteome.fasta", "fasta"))
results = SearchIO.parse("scan_results.tsv", "hmmer3-tab")

for result in results:
	for protein in proteome:
		if protein.id == result.id:
			protein.name = result.hit_keys

for protein in proteome:
	

SeqIO.write(proteome, "catproteome_scan.fasta", "fasta")
