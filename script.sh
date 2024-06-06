#!/bin/bash

# This script is the first version of NanoASV
# Authors : Arthur Cousson, Frederic Mahe
# 08/03/2023

START=$(date +%s) #Set the clock for timer

#/usr/games/cowsay -TU NanoASV is a workflow created by Arthur Cousson with useful contributions from Frederic Mahe and Enrique Ortega-Abbud. Hope this will help you analyse your data. && /usr/games/cowsay -f dragon Death To Epi2Me !
#echo "NanoASV is a workflow created by Arthur Cousson with useful contributions from Frederic Mahe and Enrique Ortega-Abbud. Don't forget to cite NanoASV and its dependencies if it helps you treating your sequencing data." #&& /usr/games/cowsay -f dragon Death To Epi2Me !

#Log system and error handling *********************************************************************************************
# LOG_FILE="NanoASV_log.txt"
# exec > >(tee -a $LOG_FILE) 2>&1
#set -e
#***************************************************************************************************************************


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
DEFAULT_TREE=1
DEFAULT_DOCKER=0
DEFAULT_R_STEP_ONLY=0
DEFAULT_METADATA=${DIR}

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
    --no-r-cleaning)
      R_CLEANING=0
      shift
      ;;
    --subsampling)
      SUBSAMPLING="$2"
      shift
      shift
      ;;
    --notree)
      TREE=0
      shift
      ;;
    --docker)
      DOCKER=1
      shift
      ;;
    --ronly)
      R_STEP_ONLY=1
      shift
      ;;
    --version)
      echo "NanoASV 1.0"
      shift
      ;;
    --metadata)
      METADATA="$2"
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
NUM_PROCESSES="${NUM_PROCESSES:-$DEFAULT_NUM_PROCESSES}"
R_CLEANING="${R_CLEANING:-$DEFAULT_R_CLEANING}"
SUBSAMPLING="${SUBSAMPLING:-$DEFAULT_SUBSAMPLING}"
TREE="${TREE:-$DEFAULT_TREE}"
DOCKER="${DOCKER:-$DEFAULT_DOCKER}"
SUBSAMPLING=$((SUBSAMPLING * 4))
R_STEP_ONLY="${R_STEP_ONLY:-$DEFAULT_R_STEP_ONLY}"
METADATA="${METADATA:-$DEFAULT_METADATA}"
#***************************************************************************************************************************


#***************************************************************************************************************************
# Check if the required binaries are correctly installed
/bin/which mafft > /dev/null || \
    { echo "mafft is not there. Please reinstall" ; exit 1 ; }

/bin/which /opt/chopper > /dev/null || \
    { echo "chopper is not there. Please reinstall" ; exit 1 ; }

/bin/which porechop > /dev/null || \
    { echo "porechop is not there. Please reinstall" ; exit 1 ; }

/bin/which bwa > /dev/null || \
    { echo "bwa is not there. Please reinstall" ; exit 1 ; }

/bin/which samtools > /dev/null || \
    { echo "samtools is not there. Please reinstall" ; exit 1 ; }

/bin/which FastTree > /dev/null || \
    { echo "FasTree is not there. Please reinstall" ; exit 1 ; }

/bin/which Rscript > /dev/null || \
    { echo "R is not there. Please reinstall" ; exit 1 ; }



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

## Create temporary directory ***********************************************************************************************
# date
# echo Creating temporary directory at /tmp/
TMP="$(mktemp --directory || exit 1)"

#****************************************************************************************************************************
if [[ "${DOCKER}" -eq 1 ]]; then
#Docker version ************************************************************************************************************
mkdir --parents \
    ${DIR}/${OUT}/Results/{ASV,Tax,Unknown_clusters,Phylogeny,Exact_affiliations,Rdata} 2> /dev/null
OUTPWD=${DIR}/${OUT}
fi

#***************************************************************************************************************************
if [[ "${DOCKER}" -eq 0 ]]; then
#Singularity version *******************************************************************************************************
mkdir --parents \
    ${OUT}/Results/{ASV,Tax,Unknown_clusters,Phylogeny,Exact_affiliations,Rdata} 2> /dev/null
OUTPWD=$(pwd)/${OUT}
fi
#***************************************************************************************************************************

#R Step Only if problem *********************************************************************************************
if [ "$R_STEP_ONLY" -eq 1 ]; then
##Production of phyloseq object *************************************************************************************
Rscript /script.r $DIR $OUTPWD $R_CLEANING $TREE

#********************************************************************************************************************
declare -i TIME=$(date +%s)-$START
#********************************************************************************************************************
echo "Data treatment is over."
echo "NanoASV Rstep took $TIME seconds to perform."
echo "Don't forget to cite NanoASV and its dependencies if it helps you treating your sequencing data."
#********************************************************************************************************************
exit
fi


## Concatenation of fastq files *********************************************************************************************

cat_files() {
  BARCODE_DIR="$1"
  # Extract the barcode from the directory name
  BARCODE=$(basename "${BARCODE_DIR}")
  # Concatenate all fastq.gz files in the barcode directory
  zcat "${DIR}/${BARCODE_DIR}"/*.fastq.gz | gzip > "${TMP}/${BARCODE}.fastq.gz"
}

export -f cat_files  # Export the function so that it can be used in parallel

echo "Step 1/9 : Concatenation"
find "${DIR}" -maxdepth 1 -type d -name "barcode*" | env TMP="${TMP}" QUAL="${QUAL}" MINL="${MINL}" MAXL="${MAXL}" ID="${ID}" \
  parallel -j "${NUM_PROCESSES}" cat_files

#***************************************************************************************************************************

## Define function to filter files******************************************************************************************

filter_file() {
  (
  filename=$(basename "$1")
  output_file="FILTERED_$filename"
  zcat "$1" | /opt/chopper -q "${QUAL}" -l "${MINL}" --maxlength "${MAXL}" 2> /dev/null | gzip > "${TMP}/${output_file}"
  )
}

export -f filter_file

#***************************************************************************************************************************


## Filtering sequences based on quality with Chopper **********************************************************************

# Iterate over the files in parallel
echo "Step 2/9 : Filtering with Chopper"
find "${TMP}" -maxdepth 1 -name "barcode*.fastq.gz" | env TMP="${TMP}" QUAL="${QUAL}" MINL="${MINL}" MAXL="${MAXL}" ID="${ID}"\
  parallel -j "${NUM_PROCESSES}" filter_file  
#echo Unfiltered files are being deleted
rm ${TMP}/barcode*.fastq.gz

# #***************************************************************************************************************************


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

echo "Step 3/9 : Adapter trimming with Porechop"
# Iterate over the files in parallel
find "${TMP}" -maxdepth 1 -name "FILTERED*.fastq.gz" | env TMP="${TMP}" QUAL="${QUAL}" MINL="${MINL}" MAXL="${MAXL}" \
ID="${ID}"  parallel -j "${NUM_PROCESSES}" chop_file  

#echo Filtered datasets are being deleted
#rm ${TMP}/NONCHIM*
rm ${TMP}/FILTERED*

# Subsampling

echo "Step 4/9 : Subsampling"


(cd ${TMP}
 for CHOPED_FILE in CHOPED*.fastq.gz ; do
     zcat "${CHOPED_FILE}" | head -n "${SUBSAMPLING}"  > "SUB_${CHOPED_FILE}"
         #echo ${CHOPED_FILE} sub-sampled
 done
)

#echo Full size datasets are being deleted
rm ${TMP}/CHOPED*
#***************************************************************************************************************************

# ## Chimera detection *******************************************************************************************************
# # Chimera detection function definition
# chimera_detection() {
#   (
#   #echo Chimera detection step
#   filename=$(basename "$1")
#   chimera_out="NONCHIM_$filename"
#   vsearch --uchime_denovo $1 --nonchimeras "${TMP}/${chimera_out}" 2> /dev/null
#   #echo ${chimera_out} chimera removed
#   )
# }
# export -f chimera_detection

echo "Step 5/9 : Chimera detection with vsearch - INACTIVATED"
# #Iterate in parallel
# find "${TMP}" -maxdepth 1 -name "FILTERED*.fastq.gz" | env TMP="${TMP}" QUAL="${QUAL}" MINL="${MINL}" MAXL="${MAXL}" ID="${ID}"\
#   parallel -j "${NUM_PROCESSES}" chimera_detection  

# #echo Filtered datasets are being deleted
# rm ${TMP}/FILTERED*

# #***************************************************************************************************************************


# Bwa alignments ***********************************************************************************************************

 SILVA="/database/SILVA_138.1_SSURef_tax_silva.fasta.gz"
 TAX="/database/Taxonomy_SILVA138.1.csv"
 DB="/database"
#***************************************************************************************************************************

# Define a function to process each file
process_file() {
    FILE="$1"
    #date
    #echo "${FILE} alignment"
    filename=$(basename "$1")
    bwa mem ${DB}/SILVA_IDX "${FILE}" 2> /dev/null > "${FILE}.sam"
    #echo samtools start
    outsamtools_file="Unmatched_$filename"
    output_file="ASV_abundance_$filename"
    samtools fastq -f 4 "${FILE}.sam" 2> /dev/null > ${TMP}/${outsamtools_file} 
    grep -v '^@' ${FILE}.sam | grep -v '[[:blank:]]2064[[:blank:]]' | grep -v '[[:blank:]]2048[[:blank:]]' | tee >(cut -f 1,2,3 > \
     "${FILE}_Exact_affiliations.tsv") | cut -f3 | sort | uniq -c | awk '$1 != 0' | sort -nr > ${TMP}/${output_file}.tsv
    sed -i 's/^[[:space:]]*//' ${TMP}/${output_file}.tsv
    grep -o '[^ ]\+$' ${TMP}/${output_file}.tsv > "${TMP}/${filename}_ASV_list.tsv"
    #echo "${FILE} taxonomy export"
    barcode_number=$(echo "$filename" | sed -E 's/.*barcode([0-9]+).*\.fastq.gz/\1/')
    output_tax="Taxonomy_barcode${barcode_number}.csv"
    grep -f "${TMP}/${filename}_ASV_list.tsv" "${TAX}" > ${TMP}/${output_tax}
}

# Export the function
export -f process_file
#***************************************************************************************************************************
echo "Step 6/9 : Reads alignements with bwa against SILVA_138.1"
# Iterate over the files in parallel
find "${TMP}" -maxdepth 1 -name "SUB_CHOPED_FILTERED_barcode*.fastq.gz" | env DB="${DB}" TMP="${TMP}" QUAL="${QUAL}" \
MINL="${MINL}" MAXL="${MAXL}" ID="${ID}" SILVA="${SILVA}" TAX="${TAX}" parallel -j "${NUM_PROCESSES}" process_file

#***************************************************************************************************************************

# Homogenization of exact affiliations file names **************************************************************************
(cd ${TMP}
for file in SUB_CHOPED_FILTERED_barcode*.fastq.gz_Exact_affiliations.tsv; do
    barcode_number=$(echo "$file" | sed -E 's/.*barcode([0-9]+).*\.tsv/\1/')
    new_file="barcode${barcode_number}_exact_affiliations.tsv"
    mv "$file" "$new_file"
done

#***************************************************************************************************************************

# Homogeneization of ASV table names ***************************************************************************************
for file in ASV_abundance_SUB_CHOPED_FILTERED_barcode*.fastq.gz.tsv; do
    barcode_number=$(echo "$file" | sed -E 's/.*barcode([0-9]+).*\.tsv/\1/')
    new_file="barcode${barcode_number}_abundance.tsv"
    mv "$file" "$new_file"
done
)
#***************************************************************************************************************************

ls ${TMP}

# Clustering step **********************************************************************************************************

# This function to hemomogeneize names
(cd ${TMP}
for file in Unmatched_SUB_CHOPED_FILTERED_barcode*.fastq.gz; do
    if [ -e "$file" ]; then
    newname=$(echo "$file" | sed 's/Unmatched_SUB_CHOPED_FILTERED_barcode\([0-9]\+\)\.fastq.gz/barcode\1_unmatched.fastq.gz/')
    mv "$file" "$newname"
    fi
done
)
#***************************************************************************************************************************

ls ${TMP}

#This function to add barcode identifier to fasta header to retrieve abundance after clustering ****************************
(cd ${TMP}
for file in barcode*_unmatched.fastq.gz; do
    if [ -e "$file" ]; then
  sample=$(echo "$file" | \
  sed 's/barcode\(.*\)_unmatched.fastq.gz/\1/');\
  awk '{if (NR%4==1) {sub("^@", "@"); print $0 ";barcodelabel=barcode'"$sample"'"} else print $0}' "$file" >\
  "$file.tmp" && mv "$file.tmp" "$file"; 
  fi
done
)
#***************************************************************************************************************************

# Vsearch Unknown sequences clustering step ********************************************************************************
UNIQ_ID=uuidgen
(cd ${TMP}
zcat barcode*_unmatched.fastq.gz > seqs 2> /dev/null
# Check if seqs is not empty
if [ -s "seqs" ]; then
echo "Step 7/9 : Unknown sequences clustering with vsearch"

vsearch \
        --cluster_size seqs \
        --id 0.7 \
        --relabel ${UNIQ_ID}_Unknown_cluster_ \
        --sizeout \
        --otutabout unknown_clusters.tsv \
        --biomout unknown_clusters.biom \
        --clusterout_id \
        --clusterout_sort \
        --consout Consensus_seq_OTU.fasta \
        #--randseed 666
rm seqs

else 
echo "Step 7/9 : Skipped - no unknown sequence"
fi
)
# Create phylogeny with MAFFT and FastTree *********************************************************************************

echo "Step 8/9 : Phylogeny with MAFFT and FastTree"

## Get every identified ASV ID

if [ "$TREE" -eq 1 ]; then
(cd ${TMP}
cat *_ASV_list.tsv | sort -u > ID_ASV
zcat ${SILVA} | grep -A 1 -f ID_ASV | grep -v "^--" > ALL_ASV.fasta
if [ -e "Consensus_seq_OTU.fasta" ]; then
cat ALL_ASV.fasta Consensus_seq_OTU.fasta > ALL_ASV_OTU.fasta
else 
cat ALL_ASV.fasta > ALL_ASV_OTU.fasta
fi

## MAFFT alignement ********************************************************************************************************
mafft --thread "${NUM_PROCESSES}" ALL_ASV_OTU.fasta > ALL_ASV.aln 2> /dev/null

## FastTree ****************************************************************************************************************
FastTree -nt ALL_ASV.aln > ASV.tree 2> /dev/null
)
fi

## Export results **********************************************************************************************************
(cd ${TMP}
mv *_abundance.tsv ${OUTPWD}/Results/ASV/
mv Taxonomy*.csv ${OUTPWD}/Results/Tax/
mv Consensus_seq_OTU.fasta unknown_clusters.tsv unknown_clusters.biom  ${OUTPWD}/Results/Unknown_clusters/ 2> /dev/null
mv *_exact_affiliations.tsv ${OUTPWD}/Results/Exact_affiliations/
mv ASV.tree ${OUTPWD}/Results/Phylogeny/
)
rm -r ${TMP}
#***************************************************************************************************************************

##Production of phyloseq object ********************************************************************************************
echo "Step 9/9 : Phylosequization with R and phyloseq"
Rscript /script.r $DIR $OUTPWD $R_CLEANING $TREE $METADATA 2> /dev/null

#***************************************************************************************************************************
declare -ir TIME=$(( $(date +%s) - ${START} ))
#***************************************************************************************************************************
echo "Data treatment is over."
echo "NanoASV took ${TIME} seconds to perform."
echo "Don't forget to cite NanoASV and its dependencies if it allows you to treat your data."
