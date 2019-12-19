# Fasta2HMMs

Bash runner and forked python script <br/>
SCRIPT:  profileHMMsFromFASTA.py <br/>
AUTHOR:  Peter Skewes-Cox <br/>
UPDATED:  February 2014 <br/>
FORKED: From "vFam", by Uri Neri Dec2019 <br/>
Fork LOG: <br/>
    A.Replaced NCBI blast with Diamond. New option [-S <bool>] sets --more-sensitive mode. <br/>
    B.Modified func. batchMuscleCall to use GNU Parl (new func. GNU_Parl_batchMuscleCall).  <br/>
    C.Enabled multi-threading when calling CD-HIT.  <br/>
    D.Added word length = 2 to the CD-HIT call (might be important if CD-HIT will be used for more than collapsing sequences). <br/>
    E.Added -P option to enable the polyprotein filtering heuristic step (*Recommend trimming the seqs to the core beforehand, so Defulat == False) <br/> 
    F.Added -M option to specify memory in Mb[default = 4800] <br/>
    
