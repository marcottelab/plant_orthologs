# plant_orthologs
Code, script, and notes for creating ortholog groups for plants (arabidopsis, rice, wheat germ, fern, etc). Requires HMMER3 and Biopython 1.66.

###Outline for generating orthogroup tables
1. Break proteome into chunks
  * Run proteome_breaker.py
  * 20,000 sequences/file is a good size
2. Hmmscan
  * This is the most expensive task, run on tacc
  * Run hmmer_scan.sh on each chunk of the proteome
    * Specify -n 1, -c 16 in header, hmmscan will automatically use all 16 cores for multithreaded parallelization
    * On stampede this scanned ~1000 sequences/hour. 
3. Process scan
  * Write a bash file running any of the process programs on each scan result from each chunk of the proteome. Add to the same output file.

###Scripts and command line examples for arabidopsis
**hmmer_scan.sh:** Takes aproteome.fasta, directory for output, name for a new hmmdatabase or path to an existing hmmdatabase, name for results file, and optional path to a directory containing hmms if an existing hmmdatabase was not provided. Returns compressed hmm database if one was not provided and hmmscan text output. Searchs each protein sequence against an hmm database.
  `bash scripts/hmmer_scan.sh input_data/arath/uniprot-proteome%3AUP000006548.0.fasta output_data/euNOG_arath/ input_data/euNOG_pressed/euNOG_hmms scan0.txt`

**process_edist.py:** Takes aproteome.fasta, hmmerscan results, level of orthologous groups searched against, name and path for output (can make new file or add to existing). Returns a space serperated text file with level, proteinID, orthogroup1 (top hit), evalue1, orthogroup2 (2nd hit), evalue2 and proteomeID (from proteome file name). If no hits orthogroup1 is listed as proteinID. For analysis of evalue distributions.
  `python scripts/process_edist.py input_data/uniprot-proteome%3AUP000006548.fasta output_data/euNOG_arath/scan0.txt eukaryotes output_data/euNOG_arath/arath_edist.txt`

**process_small.py:** Takes aproteome.fasta, hmmerscan results, level of orthologous groups searched against, name and path for output (can make new file or add to existing). Returns a space seperated text file with rank, level, proteinID, orthogroupID, evalue, and proteomeID (from proteome file name). Rank 0 orthogroup is proteinID when a protein has no significant hits (evalue>0.01). Rank 1 orthogroup is the top hit hmm if significant. Rank 2 orthogroup is the second hit, included if the first and second evalues are within 10fold difference.
  `python scripts/process_small.py input_data/uniprot-proteome%3AUP000006548.fasta output_data/euNOG_arath/scan0.txt eukaryotes output_data/euNOG_arath/arath_small.txt`

**process_tot.py:** Takes aproteome.fasta, hmmerscan results, level of orthologous groups searched against, name and path for output (can make new file or add to existing). Returns a space seperated text file with rank, level, proteinID, orthogroupID, evalue, and proteomeID (from proteome file name). Rank 0 orthogroup is proteinID when a protein has no significant hits (evalue>0.01). All hits with an e-value<=0.01 yield an entry, rank is based on order in hmmer results.
  `python scripts/process_tot.py input_data/uniprot-proteome%3AUP000006548.fasta output_data/euNOG_arath/scan0.txt eukaryotes output_data/euNOG_arath/arath_tot.txt`

**proteome_breaker.py:** Takes aproteome.fasta, number of sequences per output file, and directory to store output files. Returns a group of fasta files with specified number of sequences each(except for the last file which will likely be short).
  `python scripts/proteome_breaker.py input_data/uniprot-proteome%3AUP000006548.fasta 20000 input_data/arath/`
  
