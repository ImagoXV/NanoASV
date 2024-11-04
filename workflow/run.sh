#/bin/bash

#This script is a wrapper to launch NanoASV snakemake version with ease
conda activate NanoASV
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
        -i|--id-vsearch)
            ID="$2"
            shift
            shift
            ;;
        -ab|--minab)
            MINAB="$2"
            shift
            shift
            ;;
        -p|--num-process)
            NUM_PROCESSES="$2"
            shift
            shift
            ;;
        -db|--database)
            DATABASE="$2"
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
        -v|--version)
            echo "NanoASV 1.0 - https://github.com/ImagoXV/NanoASV - Arthur Cousson and Frederic Mahe"
            exit
            shift
            ;;
        -h|--help)
            cat /help.txt
            exit
            shift
            ;;
        --metadata)
            METADATA="$2"
            shift
            shift
            ;;
        --requirements)
            cat requirements.txt
            exit
            shift
            ;;
        *)
            echo "Unknown option: $1"
            cat /help.txt
            exit 1
            shift
            ;;
    esac
done

#***************************************************************************************************************************
# Manual entries - Arguments
# Set default values
DEFAULT_QUAL=8
DEFAULT_MINL=1300
DEFAULT_MAXL=1700
DEFAULT_ID=0.7
DEFAULT_NUM_PROCESSES=1
DEFAULT_R_CLEANING=1
DEFAULT_MINAB=5
DEFAULT_SUBSAMPLING=50000
DEFAULT_TREE=1
DEFAULT_DOCKER=0
DEFAULT_R_STEP_ONLY=0
DEFAULT_METADATA=${DIR}
#DEFAULT_DATABASE="/database/SILVA_138.1_SSURef_tax_silva.fasta.gz"
#***************************************************************************************************************************
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
SUBSAMPLING=$((SUBSAMPLING * 4))
R_STEP_ONLY="${R_STEP_ONLY:-$DEFAULT_R_STEP_ONLY}"
METADATA="${METADATA:-$DEFAULT_METADATA}"
MINAB="${MINAB:-$DEFAULT_MINAB}"

#***************************************************************************************************************************
# Check if DIR is empty and no default value is provided
if [[ -z $DIR ]]; then
    cowpy -e dead "Error: -d needs an argument, I don't know where your sequences are." >&2
    conda deactivate
    exit 1
fi
# Check if OUT is empty and no default value is provided
if [[ -z $OUT ]]; then
    cowpy -e dead "Error: -o needs an argument. You don't want me to print to stdout" >&2
    conda deactivate
    exit 1
fi

#Run the pipeline


ls $NANOASV_PATH

echo $(pwd)

snakemake -p -s "${NANOASV_PATH}"/workflow/snakefile \
    --config \
        QUAL=$QUAL \
        MINL=$MINL \
        MAXL=$MAXL \
        ID=$ID \
        NUM_PROCESSES=$NUM_PROCESSES \
        R_CLEANING=$R_CLEANING \
        MINAB=$MINAB \
        SUBSAMPLING=$SUBSAMPLING \
        TREE=$TREE \
        INPUT_DIR=$DIR \
        METADATA=$METADATA \
        DATABASE=$DATABASE

echo "OK"

conda deactivate 