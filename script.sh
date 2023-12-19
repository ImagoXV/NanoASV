#!/bin/bash

# This script is the first docker version of NanoASV
# Authors : Arthur Cousson, Frederic Mahe
# 08/03/2023

START=$(date +%s) #Set the clock for timer

#/usr/games/cowsay -TU NanoASV is a workflow created by Arthur Cousson with useful contributions from Frederic Mahe and Enrique Ortega-Abbud. Hope this will help you analyse your data. && /usr/games/cowsay -f dragon Death To Epi2Me !
#echo NanoASV is a workflow created by Arthur Cousson with useful contributions from Frederic Mahe and Enrique Ortega-Abbud. Hope this will help you analyse your data. #&& /usr/games/cowsay -f dragon Death To Epi2Me !

#***************************************************************************************************************************
# Manual entries - Arguments
# Set default values
DEFAULT_QUAL=8
DEFAULT_MINL=1300
DEFAULT_MAXL=1700
DEFAULT_ID=0.7
DEFAULT_NUM_PROCESSES=6
DEFAULT_R_CLEANING=1
DEFAULT_MINAB=0
DEFAULT_SUBSAMPLING=10000000
DEFAULT_NUM_PROCESSES=6
DEFAULT_TREE=0

#***************************************************************************************************************************
# Read the arguments passed to the script
while [[ $# -gt 0 ]]; do
  key="$1"
  case $key in
    -d|--dir)
      DIR="$2"
      shift
      shift
      ;;
    -o|--out)
      OUT="$2"
      shift
      shift
      ;;    
    -q|--quality)
      QUAL="$2"
      shift
      shift
      ;;
    -l|--minlength)
      MINL="$2"
      shift
      shift
      ;;
    -L|--maxlength)
      MAXL="$2"
      shift
      shift
      ;;
    -i|--id_vsearch)
      ID="$2"
      shift
      shift
      ;;
    -p|--num_process)
      NUM_PROCESSES="$2"
      shift
      shift
      ;;
    --r_cleaning)
      R_CLEANING="$2"
      shift
      shift
      ;;
    --subsampling)
      SUBSAMPLING="$2"
      shift
      shift
      ;;
    --tree)
      TREE="$2"
      shift
      shift
      ;;
    --version)
      echo "NanoASV 1.0"
      shift
      ;;
    
    *)
      echo "Unknown option: $1"
      shift
      ;;
  esac
done

# Assign default values if variables are empty
#DIR="/data"
QUAL="${QUAL:-$DEFAULT_QUAL}"
MINL="${MINL:-$DEFAULT_MINL}"
MAXL="${MAXL:-$DEFAULT_MAXL}"
ID="${ID:-$DEFAULT_ID}"
NUM_PROCESSES="${NUM_PROCESSES:-$DEFAULT_NUM_PROCESSES}"
R_CLEANING="${R_CLEANING:-$DEFAULT_R_CLEANING}"
SUBSAMPLING="${SUBSAMPLING:-$DEFAULT_SUBSAMPLING}"
TREE="${TREE:-$DEFAULT_TREE}"

#***************************************************************************************************************************

#***************************************************************************************************************************
# Check if DIR is empty and no default value is provided
if [[ -z $DIR ]]; then
  /usr/games/cowsay -d "Error: -d needs an argument, I don't know where your sequences are." >&2
  exit 1
fi
# Check if OUT is empty and no default value is provided
if [[ -z $OUT ]]; then
  /usr/games/cowsay -d "Error: -o needs an argument. You don't want me to print to stdout" >&2
  exit 1
fi
#***************************************************************************************************************************

## Create temporary directory ***********************************************************************************************
# date
# echo Creating temporary directory at /tmp/
mkdir /tmp/.tmp_NanoASV
TMP="/tmp/.tmp_NanoASV"


## Concatenation of fastq files *********************************************************************************************
#echo Concatenation step
(cd ${DIR}
  for BARCODE in barcode* ; do
    #  pwd
    #  ls
    #  date
    zcat -v ${BARCODE}/*.fastq.gz | gzip > ${TMP}/${BARCODE}.fastq.gz         
    #echo ${BARCODE} concatenated
 done
)
#***************************************************************************************************************************

## Define function to filter files******************************************************************************************

filter_file() {
  (
  #echo "Concerned file is $1"
  filename=$(basename "$1")
  output_file="FILTERED_$filename"
  zcat "$1" | /opt/chopper -q "${QUAL}" -l "${MINL}" --maxlength "${MAXL}" | gzip > "${TMP}/${output_file}"
  #echo "$1 filtered"
  )
}

export -f filter_file
#***************************************************************************************************************************


## Filtering sequences based on quality with NanoFilt **********************************************************************
#echo NanoFilt step
# Iterate over the files in parallel
find "${TMP}" -maxdepth 1 -name "barcode*.fastq.gz" | env TMP="${TMP}" QUAL="${QUAL}" MINL="${MINL}" MAXL="${MAXL}" ID="${ID}"\
  parallel -j "${NUM_PROCESSES}" filter_file  
#echo Unfiltered files are being deleted
rm ${TMP}/barcode*.fastq.gz
#date
#***************************************************************************************************************************

## Chimera detection *******************************************************************************************************
##vsearch --uchime_denovo FILENAME --nonchimeras FILENAME
# Chimera detection function definition
chimera_detection() {
  (
  echo Chimera detection step
  filename=$(basename "$1")
  chimera_out="NONCHIM_$filename"
  vsearch --uchime_denovo $1 --nonchimeras "${TMP}/${chimera_out}"
  echo ${chimera_out} chimera removed
  )
}
export -f chimera_detection

#Iterate in parallel
find "${TMP}" -maxdepth 1 -name "FILTERED*.fastq.gz" | env TMP="${TMP}" QUAL="${QUAL}" MINL="${MINL}" MAXL="${MAXL}" ID="${ID}"\
  parallel -j "${NUM_PROCESSES}" chimera_detection  

#echo Filtered datasets are being deleted
rm ${TMP}/FILTERED*

#***************************************************************************************************************************


## Trim adapaters with Porechop ********************************************************************************************
#echo "Porechop step"

chop_file() {
  (
  #echo "Concerned file is $1"
  filename=$(basename "$1")
  output_file="CHOPED_$filename"
  porechop --verbosity 0 -i $1 -o ${TMP}/${output_file}
  #echo "$1 choped"
  )
}

export -f chop_file
#***************************************************************************************************************************

# Iterate over the files in parallel
find "${TMP}" -maxdepth 1 -name "NONCHIM_*.fastq.gz" | env TMP="${TMP}" QUAL="${QUAL}" MINL="${MINL}" MAXL="${MAXL}" \
ID="${ID}"  parallel -j "${NUM_PROCESSES}" chop_file  

#echo Filtered datasets are being deleted
rm ${TMP}/NONCHIM*

# Subsampling
#echo Barcodes ${SUBSAMPLING} firsts quality checked sequences subsampling
#date
(cd ${TMP}
 for CHOPED_FILE in CHOPED*.fastq.gz ; do
     zcat "${CHOPED_FILE}" | head -n "${SUBSAMPLING}"  > "SUB_${CHOPED_FILE}"
         #echo ${CHOPED_FILE} sub-sampled
 done
)

#echo Full size datasets are being deleted
rm ${TMP}/CHOPED*
#***************************************************************************************************************************


# Bwa alignments ***********************************************************************************************************

 SILVA="/database/SILVA_138.1_SSURef_tax_silva.fasta.gz"
 TAX="/database/Taxonomy_SILVA138.1.csv"
 DB="/database"
#***************************************************************************************************************************

# Define a function to process each file
process_file() {
    FILE="$1"
    date
    echo "${FILE} alignment"
    filename=$(basename "$1")
    bwa mem ${DB}/SILVA_IDX "${FILE}" > "${FILE}.sam"
    echo samtools start
    outsamtools_file="Unmatched_$filename"
    output_file="ASV_abundance_$filename"
    samtools fastq -f 4 "${FILE}.sam" > ${TMP}/${outsamtools_file}
    grep -v '^@' ${FILE}.sam | grep -v '[[:blank:]]2064[[:blank:]]' | grep -v '[[:blank:]]2048[[:blank:]]' | tee >(cut -f 1,2,3 > \
     "${FILE}_Exact_affiliations.tsv") | cut -f3 | sort | uniq -c | awk '$1 != 0' | sort -nr > ${TMP}/${output_file}.tsv
    sed -i 's/^[[:space:]]*//' ${TMP}/${output_file}.tsv
    grep -o '[^ ]\+$' ${TMP}/${output_file}.tsv > "${TMP}/${filename}_ASV_list.tsv"
    echo "${FILE} taxonomy export"
    barcode_number=$(echo "$filename" | sed -E 's/.*barcode([0-9]+).*\.fastq.gz/\1/')
    output_tax="Taxonomy_barcode${barcode_number}.csv"
    grep -f "${TMP}/${filename}_ASV_list.tsv" "${TAX}" > ${TMP}/${output_tax}
}

# Export the function
export -f process_file
#***************************************************************************************************************************

# Iterate over the files in parallel
find "${TMP}" -maxdepth 1 -name "SUB_CHOPED_NONCHIM_FILTERED_barcode*.fastq.gz" | env DB="${DB}" TMP="${TMP}" QUAL="${QUAL}" \
MINL="${MINL}" MAXL="${MAXL}" ID="${ID}" SILVA="${SILVA}" TAX="${TAX}" parallel -j "${NUM_PROCESSES}" process_file

#***************************************************************************************************************************

# Homogenization of exact affiliations file names **************************************************************************
(cd ${TMP}
for file in SUB_CHOPED_NONCHIM_FILTERED_barcode*.fastq.gz_Exact_affiliations.tsv; do
    barcode_number=$(echo "$file" | sed -E 's/.*barcode([0-9]+).*\.tsv/\1/')
    new_file="barcode${barcode_number}_exact_affiliations.tsv"
    mv "$file" "$new_file"
done

#***************************************************************************************************************************

# Homogeneization of ASV table names ***************************************************************************************
for file in ASV_abundance_SUB_CHOPED_NONCHIM_FILTERED_barcode*.fastq.gz.tsv; do
    barcode_number=$(echo "$file" | sed -E 's/.*barcode([0-9]+).*\.tsv/\1/')
    new_file="barcode${barcode_number}_abundance.tsv"
    mv "$file" "$new_file"
done
)
#***************************************************************************************************************************


# Clustering step **********************************************************************************************************

# This function to hemomogeneize names
(cd ${TMP}
for file in Unmatched_SUB_CHOPED_NONCHIM__FILTERED_barcode*.fastq; do
    newname=$(echo "$file" | sed 's/Unmatched_SUB_CHOPED_NONCHIM_FILTERED_barcode\([0-9]\+\)\.fastq.gz/barcode\1_unmatched.fastq.gz/')
    mv "$file" "$newname"
done
)
#***************************************************************************************************************************

#This function to add barcode identifier to fasta header to retrieve abundance after clustering ****************************
(cd ${TMP}
for file in barcode*_unmatched.fastq; do \
sample=$(echo "$file" | \
sed 's/barcode\(.*\)_unmatched.fastq/\1/');\
awk '{if (NR%4==1) {sub("^@", "@"); print $0 ";barcodelabel=barcode'"$sample"'"} else print $0}' "$file" >\
"$file.tmp" && mv "$file.tmp" "$file"; done
)
#***************************************************************************************************************************

# Vsearch Unknown sequences clustering step ********************************************************************************
UNIQ_ID=uuidgen
(cd ${TMP}
cat barcode*_unmatched.fastq > seqs

vsearch \
        --cluster_size seqs \
        --id 0.7 \
        --relabel ${UNIQ_ID}_Unknown_cluster_ \
        --sizeout \
        --otutabout unknown_clusters.tsv \
        --biomout unknown_clusters.biom \
        --clusterout_id \
        --clusterout_sort \
        --consout Consensus_seq_OTU.fasta

rm seqs
)


# Create phylogeny with FastTree

## Get every identified ASV ID

#if [ "$TREE" -eq 1 ]; then
(cd ${TMP}
cat *_ASV_list.tsv | sort -u > ID_ASV
echo Extracting ASV SILVA fasta
zcat ${SILVA} | grep -A 1 -f ID_ASV | grep -v "^--" > ALL_ASV.fasta
cat ALL_ASV.fasta Consensus_seq_OTU.fasta > ALL_ASV.fasta

## MAFFT alignement
echo Starting MAFFT alignement
mafft --thread "${NUM_PROCESSES}" ALL_ASV.fasta > ALL_ASV.aln
echo MAFFT finished

cat ALL_ASV.aln

## FastTree
echo Starting FastTree
FastTree -nt ALL_ASV.aln > ASV.tree
echo TREE finished

cat ASV.tree
)
#fi

#***************************************************************************************************************************

# #Docker version **********************************************************************************************************
# mkdir ${DIR}/${OUT}
# mkdir ${DIR}/${OUT}/Results
# mkdir -v ${DIR}/${OUT}/Results/{ASV,Tax,Unknown_clusters,Phylogeny,Exact_affiliations,Rdata}

# OUTPWD=${DIR}/${OUT}
#***************************************************************************************************************************

#Singularity version *******************************************************************************************************
mkdir ${OUT}
mkdir ${OUT}/Results/
mkdir ${OUT}/Results/{ASV,Tax,Unknown_clusters,Phylogeny,Exact_affiliations,Rdata}

OUTPWD=$(pwd)/${OUT}
#***************************************************************************************************************************

## Export results **********************************************************************************************************
(cd ${TMP}
mv *_abundance.tsv ${OUTPWD}/Results/ASV/
mv Taxonomy*.csv ${OUTPWD}/Results/Tax/
mv Consensus_seq_OTU.fasta unknown_clusters.tsv unknown_clusters.biom  ${OUTPWD}/Results/Unknown_clusters/
mv *_exact_affiliations.tsv ${OUTPWD}/Results/Exact_affiliations/
mv ASV.tree ${OUTPWD}/Results/Phylogeny/
)
#rm -r ${TMP}
#***************************************************************************************************************************






##Production of phyloseq object ********************************************************************************************
#echo R step
Rscript /script.r $DIR $OUTPWD $R_CLEANING

#***************************************************************************************************************************
declare -i TIME=$(date +%s)-$START
#***************************************************************************************************************************
echo "Data treatment is over."
echo "NanoASV took $TIME seconds to perform."
echo "Don't forget to cite NanoASV and its dependencies if it helps you treating your sequencing data."
#***************************************************************************************************************************
