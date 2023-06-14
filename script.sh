#!/bin/bash

# This script is the first docker version of NanoASV
# Authors : Arthur Cousson, Frederic Mahe
# 08/03/2023

cowsay NanoASV is a workflow created by Arthur Cousson with useful contribution from Frederic Mahe and Enrique Ortega-Abbud. Hope this will help you analyse your data. && cowsay -f dragon Death To Epi2Me !

# Manual entries - Arguments
# Set default values
DEFAULT_QUAL=8
DEFAULT_MINL=1300
DEFAULT_MAXL=1700
DEFAULT_ID=0.7
DEFAULT_NUM_PROCESSES=6
DEFAULT_R_CLEANING=1

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
NUM_PROCESSES="${NUM_PROCESS:-$DEFAULT_NUM_PROCESSES}"
R_CLEANING="${R_CLEANING:-$DEFAULT_R_CLEANING}"


# Check if DIR is empty and no default value is provided
if [[ -z $DIR ]]; then
  echo "Error: -d needs an argument, I don't know where your sequences are."
  exit 1
fi


# Create temporary directory
######ddd#rm -vr .tmp_NanoASV

date
echo Creating temporary directory at ./.
mkdir -v .tmp_NanoASV
TMP=".tmp_NanoASV"


# Concatenation of fastq files
echo Concatenation step
(cd ${DIR} # I really need to prompt this variable as a launching option
  for BARCODE in barcode* ; do
  pwd
     date
     zcat -v ${BARCODE}/*fastq.gz > ${BARCODE}.fastq         
     echo ${BARCODE} concatenated
 done
)





# Filtering sequences based on quality with NanoFilt
echo NanoFilt step
date
cp ${DIR}/barcode*.fastq ${TMP}/
echo following stdout is ls TMP
rm ${DIR}/barcode*.fastq

chmod ugo+wrx ${TMP}/*


filter_file() {
  (
  pwd
  N_FIRST_LINES=1000000 # For optimization purposes
  echo "Concerned file is $1"
  filename=$(basename "$1")
  output_file="FILTERED_$filename"
  head -n "${N_FIRST_LINES}" "$1" | \
   NanoFilt -q "${QUAL}" -l "${MINL}" --maxlength "${MAXL}" > "${TMP}/${output_file}"
  echo "$1 filtered"
  )
}

export -f filter_file

# Iterate over the files in parallel
find "${TMP}" -maxdepth 1 -name "barcode*.fastq" | env TMP="${TMP}" QUAL="${QUAL}" MINL="${MINL}" MAXL="${MAXL}" ID="${ID}"  parallel -j "${NUM_PROCESSES}" filter_file  



echo Unfiltered files are being deleted
#rm -v ${TMP}/barcode*.fastq
date






# Chimera detection
# Work in progress

#Prochop part is working great, trying to parallel it

# Trim adapaters with Porechop


chop_file() {
  (
      echo "Concerned file is $1"
  filename=$(basename "$1")
  output_file="CHOPED_$filename"
  porechop -i $1 -o ${TMP}/${output_file} -t 16
  echo "$1 filtered"
  )
}

export -f chop_file

# Iterate over the files in parallel
find "${TMP}" -maxdepth 1 -name "FILTERED_barcode*.fastq" | env TMP="${TMP}" QUAL="${QUAL}" MINL="${MINL}" MAXL="${MAXL}" ID="${ID}"  parallel -j "${NUM_PROCESSES}" chop_file  




# Subsampling
echo Barcodes 50000 firsts quality checked sequences subsampling
date
(cd ${TMP}
 for CHOPED_FILE in CHOPED*.fastq ; do
     head -n 200000 "${CHOPED_FILE}" > "SUB_${CHOPED_FILE}"
         echo ${CHOPED_FILE} sub-sampled
 done
)


echo Full size datasets are being deleted
#rm -v ${TMP}/CHOPED*

# Bwa alignments

cp /SILVA_138.1_SSURef_tax_silva.fasta ${TMP}/SILVA.fasta

SILVA="SILVA.fasta"

chmod ugo+wrx ${TMP}/${SILVA}

# Check if the index exists
if [[ $(ls ${TMP}/*.amb 2>/dev/null | wc -l) -eq 0 ]]; then
  # Create the index
  echo Indexing SILVA
  date
  bwa index ${TMP}/${SILVA}
fi


grep ">" ${TMP}/${SILVA} | sed 's/.//' > ${TMP}/Taxonomy_SILVA138.1.csv

TAX=${TMP}/Taxonomy_SILVA138.1.csv


# Define a function to process each file
process_file() {
    FILE="$1"
    date
    echo "${FILE} alignment"
    filename=$(basename "$1")
    output_file="Abundance_ASV_$filename"
    bwa mem ${TMP}/${SILVA} "${FILE}" > "${FILE}.sam"
    samtools fastq -f 4 "${FILE}.sam" > "${FILE}_unmatched.fastq"
    grep -v '^@' "${FILE}.sam" | \
    grep -v '[[:blank:]]2064[[:blank:]]' | \ 
    grep -v '[[:blank:]]2048[[:blank:]]' | \
    tee >(cut -f 1,2,3 > "${FILE}_Exact_affiliations.tsv") | \
    cut -f3 | sort | uniq -c | awk '$1 != 1' | sort -nr > "${output_file}.tsv"

    sed -i 's/^[[:space:]]*//' "${FILE}.tsv"

    grep -o '[^ ]\+$' "${FILE}.tsv" > "${FILE}_ASV.tsv"

    echo "${FILE} taxonomy export"
    output_file="Taxonomy_$filename.csv"
    grep -f "${FILE}_ASV.txt" "${TAX}" > "${output_file}.csv"
    
}

# Export the function
export -f process_file

# Iterate over the files in parallel
find "${TMP}" -maxdepth 1 -name "SUB_CHOPEDFILTERED_barcode*.fastq.fastq" | env TMP="${TMP}" QUAL="${QUAL}" MINL="${MINL}" MAXL="${MAXL}" ID="${ID}" SILVA="${SILVA}" TAX="${TAX}" parallel -j "${NUM_PROCESSES}" process_file


# # Homogenization of exact affiliations file names
# (cd ${TMP}
# for file in SUB_CHOPEDFILTERED_barcode*.fastq.fastq_Exact_affiliations.tsv; do mv "$file" "$(echo "$file" |\
#  sed 's/SUB_CHOPEDFILTERED_\(barcode.*\)\.fastq\.fastq_Exact_affiliations\.tsv/\1_Exact_affiliations.tsv/')"; \
#  done
# )

# # Clustering step

# # This function to hemomogeneize names
# (cd ${TMP}
# for file in SUB_CHOPEDFILTERED_barcode*.fastq.fastq_unmatched.fastq; \
# do mv "$file" "$(echo "$file" | \
# sed 's/SUB_CHOPEDFILTERED_barcode\([0-9]\+\)\.fastq\.fastq_unmatched\.fastq/barcode\1_unmatched.fastq/')"; \
# done
# )
# _Taxonomy.csv
# #This function to add barcode identifier to fasta header to retrieve abundance after clustering
# (cd ${TMP}
# for file in barcode*_unmatched.fastq; do \
# sample=$(echo "$file" | \
# sed 's/barcode\(.*\)_unmatched.fastq/\1/');\
# awk '{if (NR%4==1) {sub("^@", "@"); print $0 ";barcodelabel=barcode'"$sample"'"} else print $0}' "$file" >\
# "$file.tmp" && mv "$file.tmp" "$file"; done
# )

# # Vsearch Unknown sequences clustering step

# UNIQ_ID=uuidgen
# (cd ${TMP}
# cat barcode*_unmatched.fastq > data

# vsearch \
#         --cluster_size data \
#         --minsize 8 \
#         --id 0.7 \
#         --relabel ${UNIQ_ID}_Unknown_cluster_ \
#         --sizeout \
#         --otutabout unknown_clusters.tsv \
#         --biomout unknown_clusters.biom \
#         --clusterout_id \
#         --clusterout_sort \
#         --consout Consensus_seq_OTU.fasta

# rm data

# )

# mkdir Results

# (cd ${TMP}
# mv Consensus_seq_OTU.fasta unknown_clusters.tsv unknown_clusters.biom *Exact_affiliations.tsv *_Taxonomy.csv *_ASV.tsv ~/Results/
# )

# # Production of phyloseq object
# echo R step
# Rscript -e "source('script.r')"

# # # Multialignement step

# # align_file() {
# #   FILE="$1"
# #   echo "${FILE}" Multialignement step
# #   mafft ${FILE} > ALIGN_${FILE}
# # }

# # export -f align_file

# # (cd ${TMP}

# # cat barcode*_unmatched.fastq |\
# # awk 'NR%4==1 {printf ">%s\n", substr($0,2)} NR%4==2 {print}' | \
# #  sed -e '/^$/d' > all_barcodes.fasta
 
# # cat all_barcodes.fasta Consensus_seq_OTU.fasta > data

# #  find . -maxdepth 1 -name "data" | parallel -j "${NUM_PROCESSES}" align_file {}

# #  rm data
# # )


# # # Tree construction

# # tree_construction() {
# #   FILE="$1"
# #   echo "${FILE}" 1500bp 16S phylogenetic tree construction step
# #   echo Start by homogenize fasta headers 
# #   mafft ${FILE} > ALIGN_${FILE}
# # }

# # export -f align_file

# # (cd ${TMP}
# #  find . -maxdepth 1 -name "ALIGN_*" | parallel -j "${NUM_PROCESSES}" align_file {}
# # )



# # Ongoing work 

# # The following will give abundance of clusters and their id
# # This version doesn't accept singleton nor doublon
# #grep ">" clusters.fasta | grep -v "size=1" | grep -v "size=2" | grep -o '[^>]\+$' | sed 's/.$//'
# # This version takes all of them 
# #grep ">" clusters.fasta | grep -v "size=1" | grep -v "size=2" | grep -o '[^>]\+$' | sed 's/.$//'

# # The following will homogenize the name of the files

# #for file in SUB_CHOPEDFILTERED_barcode*.fastq.fastq.tsv; do mv "$file" "$(echo "$file" | sed 's/SUB_CHOPEDFILTERED_\(barcode[0-9]*\)\.fastq\.fastq\.tsv/\1.tsv/')"; done







