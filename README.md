# Fasta2ProfileDBs 
Input is single fasta file with amino acid sequences. 
Sequences perclsutered (CD-HIT) --> "ALL vs ALL" (DIAMOND BLASTp) --> clustered (MCL) --> aligned MUSCLE --> Formatted as DBs.
Output directory should contain HMMdb (HMMER), HHMs (hh-suite) and MMseqs profile db.
*.env file specifics the given positinal argument used, which are:

   Positional arguments:
   #	Desc (suggestion)	
   1	Threads
   2	Memory in Mb
   3	output directory
   4	Input fasta file (expected *.faa)
   5	Search params used  [<E-value,score,min_alignment_coverage,qlen>] 
   6	Minimal id for preclustering sequence collapsing (0.9) 
   7	Minimal coverage for preclustering (0.75) (-AS in cd-hit, aligment coverage of the smaller seq)
   8	MCL inflation (3.6)
   9	Minimal number of sequences per cluster (2) (Do not change for now).
   10	Precluster (True)
   11	Call Diamond with Max sensitivity (True)
   12 Output clusters prefix (Iter_ID / set_ID)

    
