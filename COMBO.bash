#!/bin/bash


exit_abnormal() {
  echo "$1"
  exit "$2"
}

Help() {
  # Display Help
  echo "Usage: $0 [options] "
  echo
  echo "## Mandatory augument:"
  echo "-E               expression matrix"
  echo "-M               methylation matrix"
  #echo "-e <string>      type expression data, can be COUNT or FPKM"
  #echo "-m <string>      methylation data, can be Bvalue or Mvalue"
  echo "-i               input file with sample information"
  echo "-c               file with comparison between samples"
  echo "-h               Print this Help"
  echo "-o               output folder path"
  echo
  echo "## Optional arguments:"
  echo "# StepMiner and BooleanNet"
  echo "-d <int>         delta threshold StepMiner (default=0.5)"
  echo "-s <int>         statistic of an implication to be considered significant in BooleanNet (default=6.0)"
  echo "-P <int>         maximum p-value P of an implication to be considered significant (default=0.01)"
  echo
  echo "# Multilayer network creation"
  echo "-w <string>      metapathway, can be KEGG or KEGG_Reactome (default=KEGG)"
  echo "-f <int>         logFC threshold for expression data (default=|0.8|)"
  echo "-x <int>         logFC threshold for methylation data (default=0)"
  echo "-j <int>         adjusted pvalue threshold for methylation and expression data (default=0.05)"
  echo "-n <string>      annotation name (mutation, protein expression)"
  echo "-a               annotation table"
  echo
  echo "# Number of CPU threads"
  echo "-T <int>         threads number"
  echo
}

declare -a name
declare -a annot

while getopts "E:M:i:c:w:d::s:P:f:x:j:o:T:n:a:h" com; do
  case "${com}" in
  E) exp="${OPTARG}" ;;
  M) methyl="${OPTARG}" ;;
  i) input="${OPTARG}" ;;
  c) comparison="${OPTARG}" ;;
  w) metaph="${OPTARG}" ;;
  d) delta="${OPTARG}" ;;
  s) sig="${OPTARG}" ;;
  P) pval="${OPTARG}" ;;
  f) lgfc="${OPTARG}" ;;
  x) lgfcm="${OPTARG}" ;;
  j) adj="${OPTARG}" ;;
  o) OUTPUT="${OPTARG}" ;;
  T) THREADS="${OPTARG}" ;;
  n) name+=("${OPTARG}") ;;
  a) annot+=("${OPTARG}") ;;
  h) # display Help message
    Help
    exit
    ;;
  \?)
    exit_abnormal "Invalid option: -$OPTARG" 1
    ;;
  :)
    exit_abnormal "Option -$OPTARG requires an argument." 2
    ;;
  esac
done

#### Check parameters ####

# Check input file
[ -z "$exp" ] && exit_abnormal "Expression matrix does not exist!" 2
[ -z "$methyl" ] && exit_abnormal "Methylation matrix does not exist!" 2
[ -z "$input" ] && exit_abnormal "Sample information table does not exist!" 2
[ -z "$comparison" ] && exit_abnormal "Comparison table beetween samples does not exist!" 2
[ -z "$OUTPUT" ] && exit_abnormal "Specify Output file!" 4

# Check threads number
[ -z "$THREADS" ] && THREADS=1

# Check StepMiner and BooleanNet value and set default value
[ -z "$delta" ] && delta=0.5
[ -z "$sig" ] && sig=6
[ -z "$pval" ] && pval=0.01

# Check multilayer creation information
if [ -z "$metaph" ]; then
  metaph="KEGG"
  echo "The methapathway used for the multilayer generation is based on $metaph"
elif [ "$metaph" != "KEGG" ] || [ "$metaph" != "KEGG_Reactome" ]; then
  exit_abnormal "Specify the correct metapathway: KEGG or KEGG_Reactome!" 4
else
  echo "Expression data are in $metaph format"
fi

[ -z "$lgfc" ] && lgfc=0.8
[ -z "$lgfcm" ] && lgfcm=0
[ -z "$adj" ] && adj=0.05

#Check annotation parameters

if [ "${#name[@]}" -eq "${#annot[@]}" ]; then
  echo "Annotation information correctly loaded"
else
  exit_abnormal "Annotation information length is not correct" 4
fi

#### Generation input for multilayer ####

PATH_INPUT=$OUTPUT/input_processed
PATH_BOOL=$OUTPUT/Boolean_analysis
PATH_MULTI=$OUTPUT/multilayer
PATH_QUERY=$OUTPUT/neo4j_query
BOOL_PR=PATH/TO/BOOLEAN
scriptR=/scriptR
PATH_META=PATH/TO/METAPATHWAY
PATH_NEO4J=PATH/TO/NEO4J
NEO4J_HOME=$PATH_NEO4J/neo4j-community-4.4.7

[[ ! -d $OUTPUT ]] && mkdir "$OUTPUT"
[[ ! -d $PATH_INPUT ]] && mkdir "$PATH_INPUT"
[[ ! -d $PATH_BOOL ]] && mkdir "$PATH_BOOL"
[[ ! -d $PATH_MULTI ]] && mkdir "$PATH_MULTI"
[[ ! -d $PATH_QUERY ]] && mkdir "$PATH_QUERY"

echo "Preparing input data for multilayer creation!"

if ls "$PATH_INPUT"/log2fpkm_*_SMinput.txt >/dev/null && ls "$PATH_INPUT"/Mval_*_SMinput.txt >/dev/null; then
  echo "Input data for StepMiner analysis exist!"
else
  echo "Generate input file for Boolean Analysis"
  Rscript $scriptR/matrix_DEgene_preparation.R "$exp" "$methyl" "$PATH_INPUT" "$input" "$comparison"
fi

#StepMiner and BooleanNet
cd $BOOL_PR
cond=()
for BOOL in "$PATH_INPUT"/*_SMinput.txt; do
  step_name="${BOOL##*/}"
  s_name="${step_name%_*}"
  c_name="${s_name##*_}"
  if ls "$PATH_BOOL"/"${s_name}"_BNout.txt >/dev/null; then
    echo "Boolean implication for $s_name exists"
    cond+=("$c_name")
  else
    echo "Running StepMiner for file $s_name"
    java StepMiner -i "$BOOL" -d "$delta" -o "$PATH_BOOL"/"$s_name"_SMout.txt
    echo "Running BooleanNet for file $s_name"
    java BooleanNet -i "$PATH_BOOL"/"$s_name"_SMout.txt -s $sig -p $pval -o "$PATH_BOOL"/"$s_name"_BNout.txt
    cond+=("$c_name")
  fi
done

cd 

#Preparation custom annotation
uniq=($(printf '%s\n' "${cond[@]}" | sort -u))

if [ -z "$annot" ]; then
  annot=("NULL")
  name=("NULL")
fi

printf -v joinedn '%s,' "${name[@]}"
printf -v joineda '%s,' "${annot[@]}"

#Multilayer creation
for c in "${uniq[@]}"; do
  if ls "$PATH_MULTI"/"$c"_multilayer_nodes.csv >/dev/null; then
    echo "Multilayer data for $c exists"
  else
    echo "Multilayer creation for condition $c"
    Rscript $scriptR/multilayer_neo4j_input.R "$c" "$PATH_INPUT" "$adj" "$lgfcm" "$lgfc" "${joineda%,}" "${joinedn%,}" "$PATH_BOOL" "$metaph" "$PATH_META" "$PATH_MULTI"
  fi
done

cd $PATH_QUERY

for c in "${uniq[@]}"; do
  echo "Load Multilayer for condition $c in neo4j"
  $NEO4J_HOME/bin/neo4j start
  $NEO4J_HOME/bin/neo4j-admin set-initial-password antonioNew
  $NEO4J_HOME/bin/neo4j-admin import --force --database=neo4j --nodes="$PATH_MULTI"/"$c"_multilayer_nodes.csv --relationships="$PATH_MULTI"/"$c"_multilayer_edges.csv
  echo "Query for $c"
  python3 $PATH_NEO4J/neo4j_query_multilayer.py "$c"
  #rm -r $NEO4J_HOME/data/*
  $NEO4J_HOME/bin/neo4j stop
done

cd

echo "DONE."
