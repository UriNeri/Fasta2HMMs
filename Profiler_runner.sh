# Profier runner
# -I Sets the main inflation value to <num>. This value is the main handle for affecting cluster granularity. It is usually chosen somewhere in the range [1.2-5.0]. -I 5.0 will tend to result in fine-grained clusterings, and -I 1.2 will tend to result in very coarse grained clusterings. Your mileage will vary depending on the characteristics of your data. That is why it is a good idea to test the quality and coherency of your clusterings using clm dist and clm info. This will most likely reveal that certain values of -I are simply not right for your data. The clm dist section contains a discussion of how to use the cluster validation tools shipped with mcl (see the SEE ALSO section).
# With low values for -I, like -I 1.2, you should be prepared to use more resources in order to maintain quality of clusterings.
# echo '
# USAGE:    
#   profileHMMsFromFASTA.py -f <FASTA> [-m <#>] [-p <#>] [-c <#>] [-C <bool>] [-I <#>] [-n <#>] [-a <#>] [-M <#>] [-o <string>] [-P <bool>] [-S <bool>] [-h]

# OPTIONS:
#   -f      FASTA file sequences from which to generate clusters
#   -m      minimum sequence length [default = 1]
#   -p      minumum fraction identity for initial sequence collapsing [default = 1.0]
#   -c      minimum fraction coverage for initial sequence collapsing [default = 0.0]
#   -C      impose fraction coverage heuristics for inclusion of sequences in MSAs [default = False]
#   -I      inflation number for cluster expansion in mcl [default = 2.0]
#   -n      minimum number of sequences allowed in a MSA [default = 2]
#   -a      number of cores on which to run all processes [default = 10]
#   -M      Memort in MBs [default = 4800]
#   -o      output prefix for cluster names (default cluster)
#   -P      Polyprotein filtering [default = False]
#   -S      Run the Diamond "all-vs-all" with --more-sensitive turned on [default = False]
#   -h      print help message

# '
pyscriptfile="/path/to/profileHMMsFromFASTA_Nerid.py"

THREADS=11
input_dir="/path/to/input/dir/"
input_fasta=$input_dir"/linear_uniq_Seqs4profiles_trimmed.faa"

Attempt_ID=4

output_dir="/path/to/output/dir/Attempt_""$Attempt_ID"
mkdir $output_dir
cd $output_dir

cp $input_fasta $output_dir"/input_fasta.faa"
# pre-Profiling info
x=1.1 # Expansion parameter 
min_ali_len=80 # From the original match.

# Iter profiling params
min_cola_id=0.99
min_cola_cov=1
# min_id=1 
# min_cov=0
Memory=10000
MCL_inflation=5

# Iter profiling Booleans
coverage_heuristics="False"
polyproteins="False"
Max_sensitivity="True" 

pyscriptCMD="$($pyscriptfile -f input_fasta.faa -p $min_cola_id -c $min_cola_cov -C $coverage_heuristics -I $MCL_inflation -a $THREADS -M $Memory -o cluster -P $polyproteins -S $Max_sensitivity)"

# echo $pyscriptCMD
cd msaFiles
for i in ./cluster_*.faa #Lineraize the MSAs
do
seqkit seq -w 0 $i > temp
cat temp > $i
done
rm temp
cd ..
mv cluster.hmm Attempt_"$Attempt_ID"_profiles.hmm
hmmstat Attempt_"$Attempt_ID"_profiles.hmm > Attempt_"$Attempt_ID"_profiles_hmmstat.tsv
Nclusters=$(grep "NAME" Attempt_"$Attempt_ID"_profiles.hmm -c)
echo "Enviroment_parameters: \n
Attempt_ID = $Attempt_ID \n
coverage_heuristics = $coverage_heuristics \n
polyproteins_filtering = $polyproteins \n
Diamond_more_sensitivity = $Max_sensitivity \n
min_collapsing_id = $min_cola_id \n
min_collapsing_cov = $min_cola_cov \n
MCL_inflation = $MCL_inflation \n
Num_clusters = $Nclusters \n
Original_x = $x \n
Original_min_ali_len = $min_ali_len
###pyscriptCMD = ###$pyscriptCMD
" > Attempt_"$Attempt_ID"_params.env