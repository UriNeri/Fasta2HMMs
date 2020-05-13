Profiler_script="/path/to/Fasta2ProfileDBs.sh"
THREADS=11
Memory=4800
output_dir="/path/to/output/directory/"
input_fasta="/path/to/input/fasta.faa"
params="(0.0000000001 70 0.5 300)" #Irelevant 
min_prec_id=0.90 
min_prec_cov=0.75
MCL_inflation=2
min_nseq=2
Precluster=True
Max_sensitivity=True
Cls_Prefix=set0 # set_ID

bash $Profiler_script $THREADS $Memory $output_dir $input_fasta $params $min_prec_id $min_prec_cov $MCL_inflation $Precluster $Max_sensitivity $Cls_Prefix
