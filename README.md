# Fasta2HMMs

Bash runner and forked python script
SCRIPT:  profileHMMsFromFASTA.py
AUTHOR:  Peter Skewes-Cox
UPDATED:  February 2014
FORKED: From "vFam", by Uri Neri Dec2019
Fork LOG:
    A.Replaced NCBI blast with Diamond. New option [-S <bool>] sets --more-sensitive mode.
    B.Modified func. batchMuscleCall to use GNU Parl (new func. GNU_Parl_batchMuscleCall). 
    C.Enabled multi-threading when calling CD-HIT. 
    D.Added word length = 2 to the CD-HIT call (might be important if CD-HIT will be used for more than collapsing sequences).
    E.Added -P option to enable the polyprotein filtering heuristic step (*Recommend trimming the seqs to the core beforehand, so Defulat == False)
    F.Added -M option to specify memory in Mb[default = 4800]
    
