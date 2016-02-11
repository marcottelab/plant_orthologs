# plant_orthologs
Code, script, and notes for creating ortholog groups for plants (arabidopsis, rice, wheat germ, fern, etc). Requires HMMER3 and Biopython 1.66.

###Outline for generating orthogroup tables
1. Break proteome into chunks
  * Run proteome_breaker.py
  * 20,000 sequences/file is a good size
2. Hmmscan
  * This is the most expensive task, run on TACC
  * Based on example.sbatch write sbatch files to run hmmer_scan.sh on each chunk of the proteome
    * Specify -n 1, -c 16 in header, hmmscan will automatically use all 16 cores for multithreaded parallelization
  * Stampede has hmmer 3.0, if an hmmdatabase was pressed using hmmer 3.1 it will not recognize it.
  * If hmmer 3.1 is locally installed on Stampede be sure to change hmmer_scan.sh to include the path to the local hmmer bin.
  * On Stampede this scanned ~1000 sequences/hour. 
3. Process scan
  * Chose a process program to use depending on your needs
  * Write a bash file using run_process.sh as an example to run the process program of your choice on each scan result from each chunk of the proteome. Add to the same output file.

###Scripts and command line examples for arabidopsis
**hmmer_scan.sh:** Takes fasta file containing proteome, directory for output, name for a new hmmdatabase or path to an existing hmmdatabase, name for results file, and optional path to a directory containing hmms if an existing hmmdatabase was not provided. Returns compressed hmm database if one was not provided and hmmscan text output. Searchs each protein sequence against an hmm database.
  `bash scripts/hmmer_scan.sh input_data/arath/uniprot-proteome%3AUP000006548.0.fasta output_data/euNOG_arath/ input_data/euNOG_pressed/euNOG_hmms scan0.txt`

**process_edist.py:** Takes fasta file containing proteome, name and path of hmmerscan results, level of hmms searched against (ie eukaryotes), and name and path for output (can make new file or add to existing). Returns a space seperated text file with level, proteinID, orthogroup1 (top hit), evalue1, orthogroup2 (2nd hit), evalue2 and proteomeID (from proteome file name). If no hits orthogroup1 is listed as proteinID. For analysis of evalue distributions.
  `python scripts/process_edist.py input_data/uniprot-proteome%3AUP000006548.fasta output_data/euNOG_arath/scan0.txt eukaryotes output_data/euNOG_arath/arath_edist.txt`

**process_fusions.py:** Takes fasta file containing proteome, name and path of hmmerscan results, level of hmms searched against (ie eukaryotes), and name and path for output (can make new file or add to existing). Returns a space seperated text file with level, proteinID, orthogroupIDs, and proteome. Proteins included are potential fusion proteins that have more than one significant hit with nonoverlapping regions of alignment.
  `python scripts/process_fusions.py input_data/uniprot-proteome%3AUP000006548.fasta output_data/euNOG_arath/scan0.txt eukaryotes output_data/euNOG_arath/arath_fusions.txt`

**process_small.py:** Takes fasta file containing proteome, name and path of hmmerscan results, level of hmms searched against (ie eukaryotes), and name and path for output (can make new file or add to existing). Returns a space seperated text file with rank, level, proteinID, orthogroupID, evalue, and proteomeID (from proteome file name). Rank 0 orthogroup is proteinID when a protein has no significant hits (evalue>0.01). Rank 1 orthogroup is the top hit hmm if significant. Rank 2 orthogroup is the second hit, included if the first and second evalues are within 10fold difference.
  `python scripts/process_small.py input_data/uniprot-proteome%3AUP000006548.fasta output_data/euNOG_arath/scan0.txt eukaryotes output_data/euNOG_arath/arath_small.txt`

**process_tot.py:** Takes fasta file containing proteome, name and path of hmmerscan results, level of hmms searched against (ie eukaryotes),and name and path for output (can make new file or add to existing). Returns a space seperated text file with rank, level, proteinID, orthogroupID, evalue, and proteomeID (from proteome file name). Rank 0 orthogroup is proteinID when a protein has no significant hits (evalue>0.01). All hits with an e-value<=0.01 yield an entry, rank is based on order in hmmer results.
  `python scripts/process_tot.py input_data/uniprot-proteome%3AUP000006548.fasta output_data/euNOG_arath/scan0.txt eukaryotes output_data/euNOG_arath/arath_tot.txt`

**proteome_breaker.py:** Takes fasta file containing proteome, number of sequences per output file, and directory to store output files. Returns a group of fasta files with specified number of sequences each (except for the last file which will likely be short).
  `python scripts/proteome_breaker.py input_data/uniprot-proteome%3AUP000006548.fasta 20000 input_data/arath/`

###Note on naming conventions
Each scanned proteome has a directory for scanned and processed results, named for the set of hmms scanned against and the species of the proteome. For example, euNOG_arath refers to scanning the arabidopsis proteome against eukaryotic orthogroups from EggNOG. Within this directory scan results are named scan0.txt, the number refering to the chunk of the proteome this scan covers. Orthogroup tables are named for the species and process file used. For example, arath_tot.txt refers to the arabidopsis scan results processed by process_tot.py.  

###Outline for generating LGL images from orthogroup tables
1. Generate master file from orthogroup table (output of any of the process programs) using og_lgl_master.py
  * If want to graph ranks 1 and 2 use output of process_small.py to avoid hairball.
2. Make LGL holder directory and place master file inside.
  * Each graph should have it's own holder
3. Generate .ncol and color files using og_lgl_rank.py or og_lgl_top.py.
4. Add RBG values to color files
5. Input .ncol file into LGL  

###Scripts for generating LGL images from orthogroup tables
**og_lgl_master.py:** Takes name and path of orthogroup table, and name and path for output (can make new file or add to existing). Returns a space seperated text file with orthogroupID, proteinID, rank and species. Entries of ranks 1 and 2 are included. File is sorted by orthogroupID. Master file for creating LGL images.
  `python scripts/og_lgl_master.py output_data/euNOG_arath/arath_small.txt LGL_holders/plants_rank/plants_master.txt`

**og_lgl_rank.py:** Takes directory for input and outputs, master file (output of og_lgl_master.py), and name for .ncol. Returns .ncol file that can be fed into LGL to plot entries of ranks 1 and 2. Orthogroup will be the central node with proteins that belong to that orthogroup connected to it. Returns color_rank for coloring by rank. Rank 1 is read and rank 2 is blue, need to remove quotations (vim replace) around color before using. Returns color_file for coloring by species. Need to replace species names with RBG codes.
  `python scripts/og_lgl_rank.py LGL_holders/plants_rank/ plants_master.txt plants_rank.ncol`

**og_lgl_top.py:** Takes directory for input and outputs, master file (output of og_lgl_master.py), and name for .ncol. Returns .ncol file that can be fed into LGL to plot entries of rank 1. Returns color_file for coloring by species. Need to replace species name with RBG codes (vim replace). 
  `python scripts/og_lgl_top.py LGL_holders/plants_top/ plants_master.txt plants_top.ncol`
 
