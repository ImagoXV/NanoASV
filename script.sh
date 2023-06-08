#!/bin/bash

# This script is the first docker version of NanoASV
# Authors : Arthur Cousson, Frederic Mahe
# 08/03/2023



# Manual entries - Arguments
# Set default values
DEFAULT_DIR="default/path/to/directory"
DEFAULT_QUAL=8
DEFAULT_MINL=1400
DEFAULT_MAXL=1600
DEFAULT_ID=0.7

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
    -i|--Id_vsearch)
      ID="$2"
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
OUT="${OUT:-$DEFAULT_OUT}"
QUAL="${QUAL:-$DEFAULT_QUAL}"
MINL="${MINL:-$DEFAULT_MINL}"
MAXL="${MAXL:-$DEFAULT_MAXL}"
ID="${ID:-$DEFAULT_ID}"



# Check if DIR is empty and no default value is provided
if [[ -z $DIR ]]; then
  echo "Error: -d needs an argument, I don't know where your sequences are."
  exit 1
fi



# Clean and create the Data directory
rm -vrf "$data_dir"
mkdir -vp "$data_dir"

# Create temporary directory

rm -vr ~/.tmp_NanoASV

date
echo Creating temporary directory at ~/.
mkdir ~/.tmp_NanoASV
TMP="~/.tmp_NanoASV"


#Concatenation of fastq files
echo Concatenation step
(cd ${DIR} # I really need to prompt this variable as a launching option
  for BARCODE in barcode* ; do
     date
     cat ${BARCODE}/*fastq.gz > TMP/${BARCODE}.fastq.gz         
     echo ${BARCODE} concatenated
 done
)

# Filtering sequences based on quality with NanoFilt
echo NanoFilt step
date
(cd ${TMP} # I really need to prompt this variable as a launching option
 N_FIRST_LINES=1000000 #For optimization purposes
 for FASTQ_FILE in barcode*.fastq ; do
     head -n ${N_FIRST_LINES} "${FASTQ_FILE}" | \
         NanoFilt -q ${QUAL} -l ${MINL} --maxlength ${MAXL} > "${TMP}/FILTERED_${FASTQ_FILE}" #The length bondaries help to narrow the analysis
     echo ${FASTQ_FILE} filtered
 done
)
echo Unfiltered files are being deleted
rm -v ${TMP}/barcode*.fastq.gz
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
rm -v ${TMP}FILTERED*

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

SILVA="/data/SILVA_138.1_SSURef_tax_silva.fasta"

# Check if the index exists
if [[ $(ls /data/*.amb 2>/dev/null | wc -l) -eq 0 ]]; then
  # Create the index
  echo Indexing SILVA
  bwa index ${SILVA}
fi

echo Create SILVA Taxonomy file

grep ">" /data/SILVA_138.1_SSURef_tax_silva.fasta | sed 's/.//' > /data/Taxonomy_SILVA138.1.csv

TAX=/data/Taxonomy_SILVA138.1.csv

echo Starting alignments with bwa and SILVA138.1
(cd ${TMP}
 for FILE in SUB_*.fastq ; do
     echo ${FILE} alignment
     bwa mem ${SILVA} ${FILE} > ${FILE}.sam

     samtools fastq -f 4 ${FILE}.sam > ${FILE}_unmatched.fastq # Will extract fasta file of unlmatched sequences
     
     grep -v '^@' ${FILE}.sam > ${FILE}.sam2 # Getting rid of header
     rm -v ${FILE}.sam # Remove tmp file

     grep -v '[[:blank:]]2064[[:blank:]]' ${FILE}.sam2 > ${FILE}.sam3 # To get rid of supp alignements
     rm -v ${FILE}.sam2 # Remove tmp file

     grep -v    '[[:blank:]]2048[[:blank:]]' ${FILE}.sam3 > ${FILE}.sam # To get rid of other supp alignements
     rm -v ${FILE}.sam3 # Remove tmp file

     #Pour obtenir la table d'abondance du barcode sans les singlotons
     cat ${FILE}.sam | cut -f3 | sort | uniq -c | cut -f1 | grep -v '^      1 ' | sort -nr > ${FILE}.tsv
     sed --in-place 's/^[[:space:]]*//' ${FILE}.tsv #Removing leading spaces
     cat ${FILE}.sam | cut -f 1,2,3 > ${FILE}_Exact_affiliations.tsv #To extract exact results of alignments
          
     #Scaffold names are second field delimiter will change
     grep -o '[^ ]\+$' ${FILE}.tsv > ${FILE}_ASV.txt
     echo construction de la taxonomy de ${FILE}
     grep -f ${FILE}_ASV.txt ${TAX} > ${FILE}_Taxonomy.csv
 done
)

# Clustering step

(cd ${TMP}
 for FILE in  *_umatched.fastq ; do
 echo ${FILE} clustering with vsearch
vsearch --cluster_fast ${FILE} --centroids ${FILE}_unmatched_Clustered.fasta --id 0.7 --sizeout
)