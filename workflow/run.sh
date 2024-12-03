#/bin/bash

#This script is a wrapper to launch NanoASV snakemake version with ease
source "$(conda info --base)/etc/profile.d/conda.sh"
conda init 2> /dev/null 1> /dev/null
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
        --remove-tmp)
            TMP_FILES=0
            shift
            ;;
        -v|--version)
            echo 'NanoASV 1.1 - Conda-x-Snakemake - https://github.com/ImagoXV/NanoASV - Arthur Cousson, Frederic Mahe and Ulysse Guyet '
            exit
            shift
            ;;
        -h|--help)
            cat $NANOASV_PATH/config/help.txt
            exit
            shift
            ;;
        --metadata)
            METADATA="$2"
            shift
            shift
            ;;
        --requirements)
            cat $NANOASV_PATH/config/requirements.txt
            exit
            shift
            ;;
        --dry-run)
            DRY="n"
            DIR="$NANOASV_PATH/config/MOCK"
            OUT="$NANOASV_PATH/Mock_run_OUTPUT"
            DATABASE="$NANOASV_PATH/config/MOCK/mock_references/complete-reference.fasta"
            shift
            ;;
        --mock)
            DIR="$NANOASV_PATH/config/MOCK"
            OUT="$NANOASV_PATH/Mock_run_OUTPUT"
            DATABASE="$NANOASV_PATH/config/MOCK/mock_references/complete-reference.fasta"
            NUM_PROCESSES=2
            #TMP_FILES=0
            shift
            ;;
        --model)
            MOD="$2"
            shift
            shift
            ;;
        --rerun-incomplete)
            RERUN="--rerun-incomplete"
            shift
            ;;
        --unlock)
            UNLK="--unlock"
            shift
            ;;
        --fasttree-model)
            FASTTREE_MOD="$2"
            shift
            shift
            ;;
        --sam-qual)
            SAMQ="$2"
            shift
            shift
            ;;
        *)
            echo 'Unknown option: $1'
            cat $NANOASV_PATH/config/help.txt
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
DEFAULT_DATABASE=$NANOASV_PATH/resources/SINGLELINE_SILVA_138.2_SSURef_tax_silva.fasta.gz
DEFAULT_TMP_FILES=1
DEFAULT_MOD="map-ont"
DEFAULT_RERUN=" "
DEFAULT_FASTTREE_MOD="fastest"
DEFAULT_SAMQ=0
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
TMP_FILES="${TMP_FILES:-$DEFAULT_TMP_FILES}"
MOD="${MOD:-$DEFAULT_MOD}"
DATABASE="${DATABASE:-$DEFAULT_DATABASE}"
RERUN="${RERUN:-$DEFAULT_RERUN}"
FASTTREE_MOD="${FASTTREE_MOD:-$DEFAULT_FASTTREE_MOD}"
SAMQ="${SAMQ:-$DEFAULT_SAMQ}"

mkdir -p tmp_files 
#***************************************************************************************************************************
# Check if DIR is empty and no default value is provided
if [[ -z $DIR ]]; then
    cowpy -e dead "Error: -d needs an argument, I don't know where your sequences are." >&2
    conda deactivate
    cat $NANOASV_PATH/config/help.txt
    exit 1
fi
# Check if OUT is empty and no default value is provided
if [[ -z $OUT ]]; then
    cowpy -e dead "Error: -o needs an argument. You don't want me to print to stdout" >&2
    conda deactivate
    exit 1
fi

#Metadata sanity checks ****************************************************************************************************
(cd "${METADATA}"
 #Check if metadata.csv has been provided by the user
 [[ -s metadata.csv ]] || \
     { /usr/games/cowsay -d "Error : Please provide a metadata.csv" >&2 ; cat $NANOASV_PATH/config/requirements.txt ; exit 1 ; }

 #Check if metadata is indeed a csv and has at least 3 columns (1 rownames, two data)
 awk -F "," 'NR == 1 { exit NF > 2 ? 0 : 1}' metadata.csv || \
     { echo 'ERROR: Check metadata.csv: it does not look like a csv file. Are you sure you are using coma to separate the fields? Do you have more than two columns?' ; cat $NANOASV_PATH/config/requirements.txt ; exit 1 ; }

 #Check if metadata.csv rownames structure is correct
 awk -F "," 'NR == 1 { exit $1 == "" ? 0 : 1}' metadata.csv || \
     { echo 'ERROR: First field of first line should be empty. Please check metadata.csv file structure.' ; cat $NANOASV_PATH/config/requirements.txt ; exit 1 ; }

 #Check if metadata.csv contains enough lines
 awk 'END{ exit NR > 1 ? 0 : 1}' metadata.csv || \
     { echo 'ERROR: metadata.csv: Missing header and/or data information. Too few lines.' ; cat $NANOASV_PATH/config/requirements.txt ; exit 1 ; }


 # Check if metadata barcodes are found within DIR
 cut -f1 -d "," metadata.csv | \
     tail -n +2 | \
     while read sample_name ; do
         [[ -d ${sample_name} ]] || \
             { echo 'ERROR, ${sample_name} not found. Please check metadata.csv and barcodes directories' ; cat $NANOASV_PATH/config/requirements.txt ; exit 1 ; }
     done

 #Check if number of fields is consistent is consistent accross all number of lines
 awk -F "," '{print NF}' metadata.csv | \
     sort -u | \
     awk 'END {exit NR == 1 ? 0 : 1}' || \
     { echo 'ERROR: Check metadata.csv: not all the lines have the same number of columns' ; cat $NANOASV_PATH/config/requirements.txt ; exit 1 ; }
)

#***************************************************************************************************************************
# Check if provided minimap2 model is correct
MODS=("map-ont" "map-hifi" "map-pb" "asm5" "asm10" "asm20" "splice" "splice:hq" "ava-pb" "ava-ont")

if [[ ! " ${MODS[@]} " =~ " ${MOD} " ]]; then
    echo 'ERROR: --model invalid specification.'
    echo 'Minimap2 default model is ${DEFAULT_MOD}.'
    echo 'Available models are $(echo '${MODS[@]}')'
    echo 'WARNING, changing minimap2 alignment model will have strong repercussion on data treatment.'
    echo 'Use this option carefully. Please read minimap2 documentation.'
    exit 1
fi

#***************************************************************************************************************************
# Check if provided fasttree model is correct
F_MODS=("gtr" "gamma" "map-pb" "wag" "lg" "fastest")

if [[ ! " ${F_MODS[@]} " =~ " ${F_MODS} " ]]; then
    echo 'ERROR: --fastree-model invalid specification.'
    echo 'FastTree default model is "${DEFAULT_FASTTREE_MOD}".'
    echo 'Available models are $(echo "${F_MODS[@]}")'
    exit 1
fi

#***************************************************************************************************************************
#Ensure database exists
if [ ! -s  "${DATABASE}" ]; then
    echo "No reference file found at ${DATABASE}. If you don't want to use a personal database, ensure you already downloaded and formated SILVA in correct location."
    exit 1
fi
#***************************************************************************************************************************
#Checking reference dataset fits to requirements
echo 'Checking reference format'
#Check if sequences contain uniq ID before taxonomy
if ! grep -qE '^> ?[^ ;]+ [^ ]' "$DATABASE"; then
    cat "$NANOASV_PATH/config/requirements.txt"
        echo -e "Invalid reference file fasta header.\nNanoASV will add one for you. Keep in mind that those uniq ID will not be consistent from one run to another.\nIf you like your reference dataset this way, save it from tmp_files before deleting:'tmp_files/CORRECT_HEADER_reference.fasta'\nFor better practice, we advise you to add one yourself.\nFormatting..."

    while IFS= read -r line; do
    if [[ $line == ">"* ]]; then
        uuid=$(uuidgen)
        echo ">${uuid} ${line:1}" >> tmp_files/CORRECT_HEADER_reference.fasta
    else #Manage singleleaved and interleaved fasta file
        echo "$line" >> tmp_files/CORRECT_HEADER_reference.fasta
    fi
    done < ${DATABASE}
    $DATABASE=tmp_files/CORRECT_HEADER_reference.fasta
fi
# Ensure reference file is singleleaved fasta
if [[ ! "$(basename "$DATABASE")" == SINGLELINE_SILVA_*_SSURef_tax_silva.fasta.gz ]]; then
    LINES=$(zcat -f "${DATABASE}" | wc -l)
    SEQS=$(zgrep -c "^>" "${DATABASE}") 
    if [[ $LINES -ne $(( SEQS * 2 )) ]]; then
        echo 'Interleaved reference fasta file - Formating.'
        zcat $DATABASE | awk '/^>/ {printf("%s%s\n",(NR==1)?"":RS,$0);next;} {printf("%s",$0);} END {printf("\n");}' > tmp_files/SINGLELINE_reference.fasta
        DATABASE=tmp_files/SINGLELINE_reference.fasta
    fi
fi

#***************************************************************************************************************************
#Run the pipeline

snakemake -"${DRY}"p -s "${NANOASV_PATH}"/workflow/snakefile ${RERUN} ${UNLK}\
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
        OUT=$OUT \
        METADATA=$METADATA \
        DATABASE=$DATABASE \
        NANOASV_PATH=$NANOASV_PATH \
        MOD=$MOD \
        FASTTREE_MOD=$FASTTREE_MOD \
        SAMQ=$SAMQ

#***************************************************************************************************************************
#Remove tmp files if flag is not set
if [[ "${TMP_FILES}" -eq 0 ]]; then 
    rm -r tmp_files
fi

if [[ "${DRY}" != "n" ]]; then 
    tree $OUT
fi

conda deactivate 