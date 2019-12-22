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
    D.Added word length = 2 to the CD-HIT call <br/>
    E.Added -P option to enable the polyprotein filtering heuristic step (*Not recommend if you trim the seqs prior) <br/> 
    F.Added -M option to specify memory in Mb[default = 4800] <br/>
    G.Enabled multi-threading when calling mcl. <br/>
    
