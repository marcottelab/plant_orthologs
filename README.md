# plant_orthologs
Code and script and notes for creating ortholog groups for plants (arabidopsis, rice, wheat germ, fern, etc). Requires HMMER3 and Biopython 1.66.

###master_scripts
ortho_catsearch.sh: Takes aproteome.fasta, name and path of new directory for output, name for a new hmmdatabase or path to an existing hmmdatabase, name and path of new fasta file, and optional path to a directory containing hmms if an existing hmmdatabase was not provided. Makes new output directory and calls hmmer_search.sh and cat_search.py, returns a new directory containing the compressed hmm database if one was not provided and hmmsearch table output, as well as a fasta file with all the sequences that mapped to the same hmm concatonated. Searchs each hmm against protein sequences. 

ortho_tabscan.sh: Takes aproteome.fasta, name and path of new directory for output, name for a new hmmdatabase or path to an existing hmmdatabase, level name, name and path for new text file or an existing file to add to, and optional path to a directory containing hmms if an existing hmmdatabase was not provided. Makes new ouput directory calls hmmer_scan.sh and process_scan.py, returns a new directory containing the compressed hmmdatabase if one was not provided and the hmmscan text output, as well as a text file containg level, proteinID, GroupID and ProteomeID for each protein in aproteome.fasta. Searchs each protein sequence against an hmm database. 

###scripts
hmmer_search.sh: Takes aproteome.fasta, name and path of new directory for ouput, name for a new hmmdatabase or path to an existing hmmdatabase, and optional path to a directory containing hmms if an existing hmmdatabase was not provided. Returns a new directory containing the compressed hmm database if one was not provided and hmmsearch table output. Searchs each hmm against protein sequences. 

hmmer_scan.sh: Takes aproteome.fasta, name and path of new directory for output, name for a new hmmdatabase or path to an existing hmmdatabase, and optional path to a directory containing hmms if an existing hmmdatabase was not provided. Returns compressed hmm database if one was not provided and hmmscan table output. Searchs each protein sequence against an hmm database.

cat_search.py: Takes aproteome.fasta, path to find table output from hmmer_search.sh, and name and path of new fasta file, returns a fasta file with all the sequences that mapped to the same hmm concatonated. Only works for data from hmmsearch where the hmm is the query and the protein sequences are the result.

process_scan.py: Takes aproteome.fasta, path of to find output from hmmer_scan.sh, level and name and path of new text file or existing file to add to. Returns a text file containing level, proteinID, groupID and proteomeID for all the entries in aproteome.fasta.
