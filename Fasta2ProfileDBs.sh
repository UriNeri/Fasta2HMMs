#!/bin/bash
#hostname

if [[ $# -eq 0 ]] ; then
  echo '   
   Arguments:
   #  Desc (suggestion) 
   1  Threads
   2  Memory in Mb
   3  output directory
   4  Input fasta file (expected *.faa)
   5  Search params used  [<E-value,score,min_alignment_coverage,qlen>] 
   6  Minimal id for preclustering sequence collapsing (0.9) 
   7  Minimal coverage for preclustering (0.75) (-AS in cd-hit, aligment coverage of the smaller seq)
   8  MCL inflation (3.6)
   9  Minimal number of sequences per cluster (>=2) (Do not change for now).
   10 Precluster (True)
   11 Call Diamond with Max sensitivity (True)
   12 Output clusters prefix (Iter_ID / set_ID)
'
exit
fi
##################################################################################################
THREADS=$1
Memory=$2
output_dir=$3
input_fasta=$4
params=$5 #(0.0000000001 70 0.5 300) # [<E-value,score1/score/bits,min_alignment_coverage,qlen>]
min_prec_id=$6 #0.90 
min_prec_cov=$7 #0.75
MCL_inflation=$8 # 3.6
min_nseq=$9 # 2
Precluster=$10 # True
Max_sensitivity=$11 # True
Cls_Prefix=$12 # set_ID

ini_name=$(basename $input_fasta ".faa")
mkdir $output_dir
cd $output_dir
mkdir msaFiles singletons fastaFiles HHMfiles logs HMMfiles MMseqs2_profiles Clusteringfiles NotEnoughSeqs
if [ "$Precluster" == "True" ]; 
then
  mkdir preclusters
  cd-hit -n 3 -M $Memory -T $THREADS  -i $input_fasta -o preclusters/"$ini_name"_c"$min_prec_id" -AS $min_prec_cov -c $min_prec_id 
  input_fasta="$output_dir"/preclusters/"$ini_name"_c"$min_prec_id" 
fi

cd Clusteringfiles
if [ "$Max_sensitivity" == "True" ]; 
then
  diamond blastp -q $input_fasta -o "$ini_name"_allByAll_blastp.br --db $input_fasta --outfmt 6  --more-sensitive  --threads $THREADS 
fi
if [ "$Max_sensitivity" == "False" ]; # Keep different IFs, (might change to NCBI BLASTp option later).
then
  diamond blastp -q $input_fasta -o "$ini_name"_allByAll_blastp.br --db $input_fasta --outfmt 6   --threads $THREADS 
fi
awk '{print $1,$2,$11}' "$ini_name"_allByAll_blastp.br  >  "$ini_name"_allByAll_blastp.abc
mcxload -abc "$ini_name"_allByAll_blastp.abc --stream-mirror --stream-neg-log10 -stream-tf "ceil(200)" -o "$ini_name"_allByAll_blastp.mci -write-tab "$ini_name"_allByAll_blastp.tab
mcl "$ini_name"_allByAll_blastp.mci -use-tab "$ini_name"_allByAll_blastp.tab -I $MCL_inflation -o "$ini_name"_allByAll_blastp.mcl -t $THREADS 

if [ $(which seqkit 2>/dev/null) ]; then
  echo "seqkit found, using that"
i=0
while IFS="" read -r p || [ -n "$p" ]
do
   i=$((i+1))
  cls_mems=$(printf '%s' "$p")
  cls_mems1=$(sed  's/ /'',''/g' <<< $(echo $cls_mems))
  seqkit grep $input_fasta -p "$cls_mems1" -w 0 > "$Cls_Prefix"."$i".faa
done < "$ini_name"_allByAll_blastp.mcl
else
echo "seqkit not found, using awk"
i=0
while IFS="" read -r p || [ -n "$p" ]
do
   i=$((i+1))
  cls_mems=$(printf '%s' "$p")
  sed  's/ /''\n>''/g' <<< $(echo $cls_mems) >cls_mems1
  awk 'BEGIN{RS=">";FS="\n"}NR==FNR{a[$1]++}NR>FNR{if ($1 in a && $0!="") printf ">%s",$0}' cls_mems1 $input_fasta > "$Cls_Prefix"."$i".faa
done < "$informati_name"_allByAll_blastp.mcl
rm cls_mems1
fi

singlist=$(grep ">" -F ./*.faa -c | grep ":1$"  |  sed  's/':1'/ /g' | sed  's|/\n||g')
mv $singlist ../singletons/
cat ../singletons/* > ../singletons.faa

### filter faa's with #nseq<min_nseq ###
if (( $min_nseq > 2)); # 
then
grep ">" -F ./*.faa -c |  sed  's/':'/\t/g' > nseqtbl
echo "$(echo $(awk  -v nseq=$min_nseq '$2 < nseq || NR==1' nseqtbl ))"  | sed  's/ .[0-9]* / /g' > Nseqtbl
mv $(cat Nseqtbl) ../NotEnoughSeqs/
mv ../NotEnoughSeqs/$clssn 
rm nseqtbl
fi

parallel -j"$THREADS" muscle -in {} -out {.}.msa.faa -log {.}.faa.log -quiet  ::: ./*.faa
mv ./*.log ../logs
mv ./*.msa.faa ../msaFiles/ 

cd ../
for i in ./msaFiles/*.faa
do
file_with_suffix=$(basename $i) 
file_name=$(basename $file_with_suffix ".msa.faa")
hmmbuild --informat afa -n $file_name -o logs/"$file_name".hmm.log HMMfiles/"$file_name".hmm  $i
hhmake -v 2 -name $file_name -i $i -o ./HHMfiles/"$file_name".hhm
done 
ffindex_build "$Cls_Prefix"_p.hhm "$Cls_Prefix"_p.hhm.index ./HHMfiles/*.hhm
mmseqs convertprofiledb "$Cls_Prefix"_p.hhm ./MMseqs2_profiles/"$Cls_Prefix"_MM
cat ./HMMfiles/*.hmm > "$Cls_Prefix"_p.hmm
hmmstat "$Cls_Prefix"_p.hmm > "$Cls_Prefix"_p_hmmstat.tsv
mv ./*.log ./logs/
cd msaFiles
for i in ./*.faa #Lineraize the MSAs
do
awk 'BEGIN{RS=">";FS="\n"}NR>1{seq="";for (i=2;i<=NF;i++) seq=seq""$i; print ">"$1"\n"seq}' $i > temp.faa
cat temp.faa > $i
done
rm temp.faa
cd ../

Nsingletons=$(grep ">" singletons.faa -c)
Nprofiles=$(grep "NAME" "$Cls_Prefix"_p.hmm -c)
echo "Enviroment_parameters: 
Output clusters prefix = $12 
Diamond_more_sensitivity = $Max_sensitivity 
Preclustering id >= $min_prec_id 
Preclustering coverage >= $min_prec_cov 
MCL_inflation = $MCL_inflation 
min_nseq = $min_nseq 
Nprofiles = $Nprofiles 
Nsingletons = $Nsingletons 
Original_params $(echo ${params[*]}) [<E-value,score1,min_alignment_coverage,qlen>] 
" > "$Cls_Prefix"_profiler.env

echo "Done: $Cls_Prefix"
