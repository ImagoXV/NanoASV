#!/bin/bash

# This script is the first docker version of NanoASV
# Authors : Arthur Cousson, Frederic Mahe
# 08/03/2023



# Manual entries - Arguments
# Set default values
DEFAULT_QUAL=8
DEFAULT_MINL=1400
DEFAULT_MAXL=1600
DEFAULT_ID=0.7
DEFAULT_NUM_PROCESSES=6

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
    *)
      echo "Unknown option: $1"
      shift
      ;;
  esac
done

# Assign default values if variables are empty
DIR="/data"
QUAL="${QUAL:-$DEFAULT_QUAL}"
MINL="${MINL:-$DEFAULT_MINL}"
MAXL="${MAXL:-$DEFAULT_MAXL}"
ID="${ID:-$DEFAULT_ID}"
NUM_PROCESSES="${NUM_PROCESS:-$DEFAULT_NUM_PROCESSES}"



# Check if DIR is empty and no default value is provided
if [[ -z $DIR ]]; then
  echo "Error: -d needs an argument, I don't know where your sequences are."
  exit 1
fi


# Create temporary directory
ls 
#rm -vr /data/.tmp_NanoASV

date
echo Creating temporary directory at ./.
mkdir -v .tmp_NanoASV
TMP=".tmp_NanoASV"


#Concatenation of fastq files
# echo Concatenation step
# (cd ${DIR} # I really need to prompt this variable as a launching option
#   for BARCODE in barcode* ; do
#      date
#      zcat ${BARCODE}/*fastq.gz > TMP/${BARCODE}.fastq         
#      echo ${BARCODE} concatenated
#  done
# )

# Filtering sequences based on quality with NanoFilt
echo NanoFilt step
date
cp ${DIR}/barcode*.fastq ${TMP}
echo following stdout is ls TMP
ls ${TMP}

#NanoFilt --help

(cd ${TMP} 
 pwd
 N_FIRST_LINES=1000000 #For optimization purposes
 for FASTQ_FILE in barcode*.fastq ; do
 echo Concerned file is ${FASTQ_FILE}
     head -n ${N_FIRST_LINES} ${FASTQ_FILE} | \
         NanoFilt -q ${QUAL} -l ${MINL} --maxlength ${MAXL} > FILTERED_${FASTQ_FILE} #The length bondaries help to narrow the analysis
     echo ${FASTQ_FILE} filtered
 done
)
echo Unfiltered files are being deleted
rm -v ${TMP}/barcode*.fastq
date


# Chimera detection
# Work in progress

# Trim adapaters with Porechop
echo Porechop step
date
(cd ${TMP}
 for TRIMMED_FILE in FILTERED*.fastq ; do
     porechop -i "${TRIMMED_FILE}" -o CHOPED${TRIMMED_FILE}.fastq -t 16
         echo ${TRIMMED_FILE} choped
 done
)

echo Unchoped files are being deleted
rm -v ${TMP}/FILTERED*

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
rm -v ${TMP}/CHOPED*

# Bwa alignments

SILVA="/SILVA_138.1_SSURef_tax_silva.fasta"

# Check if the index exists
if [[ $(ls *.amb 2>/dev/null | wc -l) -eq 0 ]]; then
  # Create the index
  echo Indexing SILVA
  #bwa index ${SILVA}
fi


grep ">" /SILVA_138.1_SSURef_tax_silva.fasta | sed 's/.//' > /Taxonomy_SILVA138.1.csv

TAX=/Taxonomy_SILVA138.1.csv


# Define a function to process each file
process_file() {
    FILE="$1"
    date
    echo "${FILE} alignment"
    cd "${TMP}"
    bwa mem "${SILVA}" "${FILE}" > "${FILE}.sam"
    samtools fastq -f 4 "${FILE}.sam" > "${FILE}_unmatched.fastq"
    grep -v '^@' "${FILE}.sam" | \
    grep -v '[[:blank:]]2064[[:blank:]]' | \ 
    grep -v '[[:blank:]]2048[[:blank:]]' | \
    tee >(cut -f 1,2,3 > "${FILE}_Exact_affiliations.tsv") | \
    cut -f3 | sort | uniq -c | awk '$1 != 1' | sort -nr > "${FILE}.tsv"

    sed -i 's/^[[:space:]]*//' "${FILE}.tsv"

    grep -o '[^ ]\+$' "${FILE}.tsv" > "${FILE}_ASV.txt"

    echo "Construction de la taxonomy de ${FILE}"
    grep -f "${FILE}_ASV.txt" "${TAX}" > "${FILE}_Taxonomy.csv"

}



# Export the function
export -f process_file

# Iterate over the files in parallel
find "${TMP_DIR}" -maxdepth 1 -name "SUB_*.fastq" | parallel -j "${NUM_PROCESSES}" process_file


# Clustering step

(cd ${TMP}
 find . -maxdepth 1 -name "*_unmatched.fastq" | parallel -j "${NUM_PROCESSES}" cluster_file {}

cluster_file() {
  FILE="$1"
  echo "${FILE}" clustering with vsearch
  vsearch --cluster_fast "${FILE}" --centroids "${FILE}_unmatched_Clustered.fasta" --id 0.7 --sizeout
}
