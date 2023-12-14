#!/bin/bash

# This script is the first docker version of NanoASV
# Authors : Arthur Cousson, Frederic Mahe
# 08/03/2023

START=$(date +%s) #Nombre de secondes depuis le debut d'Unix

#/usr/games/cowsay -TU NanoASV is a workflow created by Arthur Cousson with useful contributions from Frederic Mahe and Enrique Ortega-Abbud. Hope this will help you analyse your data. && /usr/games/cowsay -f dragon Death To Epi2Me !
#echo NanoASV is a workflow created by Arthur Cousson with useful contributions from Frederic Mahe and Enrique Ortega-Abbud. Hope this will help you analyse your data. #&& /usr/games/cowsay -f dragon Death To Epi2Me !

# Manual entries - Arguments
# Set default values
DEFAULT_QUAL=8
DEFAULT_MINL=1300
DEFAULT_MAXL=1700
DEFAULT_ID=0.7
DEFAULT_NUM_PROCESSES=6
DEFAULT_R_CLEANING=1
DEFAULT_MINAB=0

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
# OUT="${OUT:-/dev/stdout}"


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


# Create temporary directory
# r m -vr .tmp_NanoASV #To remove the tmp previous tmp directory 

date
echo Creating temporary directory at /tmp/
mkdir -v /tmp/.tmp_NanoASV
TMP="/tmp/.tmp_NanoASV"


# Concatenation of fastq files
echo Concatenation step
(cd ${DIR} # I really need to prompt this variable as a launching option
  for BARCODE in barcode* ; do
     pwd
     ls
     date
    zcat -v ${BARCODE}/*.fastq.gz | gzip > ${TMP}/${BARCODE}.fastq.gz         
    echo ${BARCODE} concatenated
 done
)





filter_file() {
  (
  N_FIRST_LINES=1000000 # For optimization purposes
  echo "Concerned file is $1"
  filename=$(basename "$1")
  output_file="FILTERED_$filename"
  zcat "$1" | head -n "${N_FIRST_LINES}" | \
   NanoFilt -q "${QUAL}" -l "${MINL}" --maxlength "${MAXL}" | gzip > "${TMP}/${output_file}"
  echo "$1 filtered"
  )
}

export -f filter_file

# Filtering sequences based on quality with NanoFilt
echo NanoFilt step
# Iterate over the files in parallel
find "${TMP}" -maxdepth 1 -name "barcode*.fastq.gz" | env TMP="${TMP}" QUAL="${QUAL}" MINL="${MINL}" MAXL="${MAXL}" ID="${ID}"  parallel -j "${NUM_PROCESSES}" filter_file  

ls ${TMP}

echo Unfiltered files are being deleted
rm -v ${TMP}/barcode*.fastq.gz
date


# Chimera detection
echo "Chimera detection - WORK IN PROGRESS"
# Work in progress

# # Trim adapaters with Porechop

echo "Porechop step ************************************************************************"

chop_file() {
  (
      echo "Concerned file is $1"
  filename=$(basename "$1")
  output_file="CHOPED_$filename"
  porechop -i $1 -o ${TMP}/${output_file} -t 4
  echo "$1 choped"
  )
}

export -f chop_file

# Iterate over the files in parallel
find "${TMP}" -maxdepth 1 -name "FILTERED_barcode*.fastq.gz" | env TMP="${TMP}" QUAL="${QUAL}" MINL="${MINL}" MAXL="${MAXL}" ID="${ID}"  parallel -j "${NUM_PROCESSES}" chop_file  

echo Filtered datasets are being deleted
rm -v ${TMP}/FILTERED*

# Subsampling
echo Barcodes 50000 firsts quality checked sequences subsampling
date
(cd ${TMP}
 for CHOPED_FILE in CHOPED*.fastq.gz ; do
     zcat "${CHOPED_FILE}" | head -n 200000  > "SUB_${CHOPED_FILE}"
         echo ${CHOPED_FILE} sub-sampled
 done
)

echo Full size datasets are being deleted
rm -v ${TMP}/CHOPED*

ls ${TMP}

# Bwa alignments

 SILVA="database/SILVA_138.1_SSURef_tax_silva.fasta.gz"


# # Check if the index exists
# echo
# if [[ $(ls ${TMP}/*.amb 2>/dev/null | wc -l) -eq 0 ]]; then
#   # Create the index
#   echo Indexing SILVA
#   date
#   bwa index ${TMP}/${SILVA}
#   grep ">" ${TMP}/${SILVA} | sed 's/.//' > ${TMP}/Taxonomy_SILVA138.1.csv

# fi



TAX=database/Taxonomy_SILVA138.1.csv

DB="/database"


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
    grep -v '^@' ${FILE}.sam | grep -v '[[:blank:]]2064[[:blank:]]' | grep -v '[[:blank:]]2048[[:blank:]]' | tee >(cut -f 1,2,3 > "${FILE}_Exact_affiliations.tsv") | cut -f3 | sort | uniq -c | awk '$1 != 0' | sort -nr > ${TMP}/${output_file}.tsv

    sed -i 's/^[[:space:]]*//' ${TMP}/${output_file}.tsv

    grep -o '[^ ]\+$' ${TMP}/${output_file}.tsv > "${TMP}/${filename}_ASV_list.tsv"

    echo "${FILE} taxonomy export"
    barcode_number=$(echo "$filename" | sed -E 's/.*barcode([0-9]+).*\.fastq.gz/\1/')
    output_tax="Taxonomy_barcode${barcode_number}.csv"
    grep -f "${TMP}/${filename}_ASV_list.tsv" "${TAX}" > ${TMP}/${output_tax}
    
}

# Export the function
export -f process_file

# Iterate over the files in parallel
find "${TMP}" -maxdepth 1 -name "SUB_CHOPED_FILTERED_barcode*.fastq.gz" | env DB="${DB}" TMP="${TMP}" QUAL="${QUAL}" MINL="${MINL}" MAXL="${MAXL}" ID="${ID}" SILVA="${SILVA}" TAX="${TAX}" parallel -j "${NUM_PROCESSES}" process_file


# Homogenization of exact affiliations file names
(cd ${TMP}
for file in SUB_CHOPED_FILTERED_barcode*.fastq.gz_Exact_affiliations.tsv; do
    barcode_number=$(echo "$file" | sed -E 's/.*barcode([0-9]+).*\.tsv/\1/')
    new_file="barcode${barcode_number}_exact_affiliations.tsv"
    mv "$file" "$new_file"
done



# Homogeneization of ASV table names
for file in ASV_abundance_SUB_CHOPED_FILTERED_barcode*.fastq.gz.tsv; do
    barcode_number=$(echo "$file" | sed -E 's/.*barcode([0-9]+).*\.tsv/\1/')
    new_file="barcode${barcode_number}_abundance.tsv"
    mv "$file" "$new_file"
done
)



# Clustering step

# This function to hemomogeneize names
(cd ${TMP}
for file in Unmatched_SUB_CHOPED_FILTERED_barcode*.fastq; do
    newname=$(echo "$file" | sed 's/Unmatched_CHOPED_FILTERED_barcode\([0-9]\+\)\.fastq.gz/barcode\1_unmatched.fastq.gz/')
    mv "$file" "$newname"
done

)


#This function to add barcode identifier to fasta header to retrieve abundance after clustering
(cd ${TMP}
for file in barcode*_unmatched.fastq; do \
sample=$(echo "$file" | \
sed 's/barcode\(.*\)_unmatched.fastq/\1/');\
awk '{if (NR%4==1) {sub("^@", "@"); print $0 ";barcodelabel=barcode'"$sample"'"} else print $0}' "$file" >\
"$file.tmp" && mv "$file.tmp" "$file"; done
)

# Vsearch Unknown sequences clustering step

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


# mkdir ${DIR}/${OUT}
# mkdir ${DIR}/${OUT}/Results
# mkdir ${DIR}/${OUT}/Results/Tax
# mkdir ${DIR}/${OUT}/Results/ASV
# mkdir ${DIR}/${OUT}/Results/Unknown_clusters
# mkdir ${DIR}/${OUT}/Results/Exact_affiliations
# mkdir ${DIR}/${OUT}/Results/Rdata

mkdir ${OUT}
mkdir ${OUT}/Results
mkdir ${OUT}/Results/Tax
mkdir ${OUT}/Results/ASV
mkdir ${OUT}/Results/Unknown_clusters
mkdir ${OUT}/Results/Exact_affiliations
mkdir ${OUT}/Results/Rdata

ls ${TMP}

(cd ${TMP}
mv *_abundance.tsv ${OUT}/Results/ASV/
mv Taxonomy*.csv ${OUT}/Results/Tax/
mv Consensus_seq_OTU.fasta unknown_clusters.tsv unknown_clusters.biom  ${OUT}/Results/Unknown_clusters/
mv *_exact_affiliations.tsv ${OUT}/Results/Exact_affiliations/
)


# cp -r ${DIR}/${OUT}/* ${DIR}/${OUT}/Results/Rdata

# Production of phyloseq object
echo R step
Rscript -e 'source("/script.r")'


declare -i TIME=$(date +%s)-$START

echo "Data treatment is over.********************************************************"
echo "It took $TIME seconds to perform."

echo "Don't forget to cite NanoASV if it helps you treating your sequencing data."

echo "Don't forget to cite NanoASV dependencies as well !****************************"

#rm -rf ${TMP}


































































# # Production of phyloseq object
# echo R step
# Rscript -e "source('script.r')"


# r m -r ${TMP}


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
 
# # cat all_barcodes.fasta Consensus_seq_OTU.fasta > dataSUB_CHOPEDFILTERED_barcode37.fastq.fastq.sam

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






