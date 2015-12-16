# plant_orthologs
Code and script and notes for creating ortholog groups for plants (arabidopsis, rice, wheat germ, fern, etc).

###master_scripts
ortho_catsearch.sh: Takes a directory of hmm files, aproteome.fasta, name and path of new directory for output, name and path of new fasta file, calls hmmer_search.sh and cat_search.py, returns a new directory containing the compressed hmm database and hmmsearch table output and fasta file with all the sequences that mapped to the same hmm concatonated. Searchs each hmm against protein sequences. 

ortho_tabscan.sh: 

###scripts
hmmer_search.sh: Takes a directory of hmm files, aproteome.fasta, name and path of new directory for ouput,
returns a new directory containing the compressed hmm database and hmmsearch table output. Searchs each hmm against protein sequences. 

hmmer_scan.sh: Returns compressed hmm database and hmmscan table output. File names hard coded right now, will change. Searchs each protein sequence against hmm database.

cat_search.py: Takes aproteome.fasta, path of where to find table output from hmmer_search.sh, and name and path of new fasta file, returns a fasta file with all the sequences that mapped to the same hmm concatonated. Only works for data from hmmsearch where the hmm is the query and the protein sequences are the result.

process_scan.py:
