#!/bin/bash

#This script is a wrapper to launch NanoASV snakemake version with ease
trap 'echo "ERROR: Something went wrong - Carefully read the above message to understand how to fix it. If you think this is a bug, please reach us out on github" >&2' ERR
set -e
# functions inherit the ERR trap
set -E

readonly MOCK_DIR="${NANOASV_PATH}/config/MOCK"
readonly MOCK_OUT="${NANOASV_PATH}/Mock_run_OUTPUT"
readonly MOCK_DATABASE="${MOCK_DIR}/mock_references/complete-reference.fasta"
readonly MOCK_NUM_PROCESSES=2

readonly DEFAULT_QUAL=8
readonly DEFAULT_MINL=1300
readonly DEFAULT_MAXL=1700
readonly DEFAULT_ID=0.7
readonly DEFAULT_NUM_PROCESSES=1
readonly DEFAULT_R_CLEANING=1
readonly DEFAULT_MINAB=5
readonly DEFAULT_SUBSAMPLING=50000
readonly DEFAULT_TREE=1
readonly DEFAULT_DATABASE="${NANOASV_PATH}/resources/SINGLELINE_SILVA_138.2_SSURef_tax_silva.fasta.gz"
readonly DEFAULT_TMP_FILES=1
readonly DEFAULT_MOD="map-ont"
readonly DEFAULT_FASTTREE_MOD="fastest"
readonly DEFAULT_SAMQ=0
readonly LINES_PER_FASTQ_ENTRY=4
readonly -a MINIMAP2_MODELS=("map-ont" "map-hifi" "map-pb" "asm5" "asm10" "asm20" "splice" "splice:hq" "ava-pb" "ava-ont")
readonly -a FASTTREE_MODELS=("gtr" "gamma" "wag" "lg" "fastest")

# usage: is_one_of VALUE ITEM...
is_one_of () {
    local -r value="${1}"
    shift
    local item
    for item in "$@" ; do
        [[ "${item}" == "${value}" ]] && return 0
    done
    return 1
}

check_minimap2_model () {
    is_one_of "${MOD}" "${MINIMAP2_MODELS[@]}" && return 0
    {
        echo 'ERROR: --model invalid specification.'
        echo "Minimap2 default model is ${DEFAULT_MOD}."
        echo "Available models are ${MINIMAP2_MODELS[*]}"
        echo 'WARNING, changing minimap2 alignment model will have strong repercussion on data treatment.'
        echo 'Use this option carefully. Please read minimap2 documentation.'
    } >&2
    exit 1
}

check_fasttree_model () {
    is_one_of "${FASTTREE_MOD}" "${FASTTREE_MODELS[@]}" && return 0
    {
        echo 'ERROR: --fasttree-model invalid specification.'
        echo "FastTree default model is \"${DEFAULT_FASTTREE_MOD}\"."
        echo "Available models are ${FASTTREE_MODELS[*]}"
    } >&2
    exit 1
}

check_reference_exists () {
    [[ -s "${DATABASE}" ]] && return 0
    echo "No reference file found at ${DATABASE}. If you don't want to use a personal database, ensure you already downloaded and formated SILVA in correct location." >&2
    exit 1
}

check_reference_headers () {
    # minimap2 gives an empty name to a sequence whose header starts with '> '
    zgrep -q '^> ' "${DATABASE}" || return 0
    {
        cat "${NANOASV_PATH}/config/requirements.txt"
        echo "ERROR: reference fasta headers must not have a space after '>'. First one found:"
        zcat -f "${DATABASE}" | grep -m 1 '^> '
    } >&2
    exit 1
}

# usage: add_reference_ids (updates DATABASE)
# headers without a unique ID before the taxonomy get a random uuid
add_reference_ids () {
    zgrep -qE '^>[^ ;]+ [^ ]' "${DATABASE}" && return 0
    cat "${NANOASV_PATH}/config/requirements.txt"
        echo -e "Invalid reference file fasta header.\nNanoASV will add one for you.\n\
        Keep in mind that those uniq ID will not be consistent from one run to another.\n\
        If you like your reference dataset this way, save it from tmp_files before deleting:'tmp_files/CORRECT_HEADER_reference.fasta'\n\
        For better practice, we advise you to add one yourself.\n\
        Formatting will take a while ..."

    # random uuid v4 in awk: one uuidgen call per header takes ~20 min on SILVA
    zcat -f "${DATABASE}" | \
        awk 'function hex(n,   s) {
                 s = ""
                 while (n-- > 0) s = s sprintf("%x", int(rand() * 16))
                 return s
             }
             BEGIN { srand() }
             /^>/ {
                 printf ">%s-%s-4%s-%x%s-%s %s\n", hex(8), hex(4), hex(3),
                        8 + int(rand() * 4), hex(3), hex(12), substr($0, 2)
                 next
             }
             { print } #Manage singleleaved and interleaved fasta file
            ' > tmp_files/CORRECT_HEADER_reference.fasta
    DATABASE="tmp_files/CORRECT_HEADER_reference.fasta"
}

# usage: linearize_reference (updates DATABASE)
# one sequence line per entry; the SILVA file made as in the README
# already is
linearize_reference () {
    [[ "$(basename "${DATABASE}")" == SINGLELINE_SILVA_*_SSURef_tax_silva.fasta.gz ]] && return 0
    local lines
    lines="$(zcat -f "${DATABASE}" | wc -l)"
    local sequences
    sequences="$(zgrep -c "^>" "${DATABASE}")"
    (( lines == sequences * 2 )) && return 0
    echo 'Interleaved reference fasta file - Formating.'
    zcat -f "${DATABASE}" | awk '/^>/ {printf("%s%s\n",(NR==1)?"":RS,$0);next;} {printf("%s",$0);} END {printf("\n");}' > tmp_files/SINGLELINE_reference.fasta
    DATABASE="tmp_files/SINGLELINE_reference.fasta"
}

# usage: gzip_fastq_in_barcode_dirs (updates DIR and METADATA)
# barcode directories holding uncompressed fastq files are copied to
# tmp_files/tmp_seq_dir/, with every fastq file gzipped
gzip_fastq_in_barcode_dirs () {
    find "${DIR}" -mindepth 2 -name "*.fastq" | grep -q . || return 0
    echo "NanoASV found flat fastq files. Compressing temporary files for compatibility"
    local -r tmp_seq_dir="tmp_files/tmp_seq_dir"
    mkdir -p "${tmp_seq_dir}"
    local barcode_dir
    local tmp_barcode_dir
    local file
    for barcode_dir in "${DIR}"/barcode*/; do
        [[ -d "${barcode_dir}" ]] || continue
        tmp_barcode_dir="${tmp_seq_dir}/$(basename "${barcode_dir}")"
        mkdir -p "${tmp_barcode_dir}"
        echo "${tmp_barcode_dir}"
        for file in "${barcode_dir}"/*; do
            if [[ -f "${file}" ]]; then
                if [[ "${file}" == *.fastq.gz ]]; then
                    ln -s "$(realpath "${file}")" "${tmp_barcode_dir}/"
                elif [[ "${file}" == *.fastq ]]; then
                    gzip -c "${file}" > "${tmp_barcode_dir}/$(basename "${file}").gz"
                fi
            fi
        done
    done
    cp "${METADATA}/metadata.csv" "${tmp_seq_dir}/metadata.csv"
    DIR="${tmp_seq_dir}"
    METADATA="${DIR}"
}

# usage: prepare_input_layout (updates DIR and METADATA)
# expects barcode directories; fastq files given without barcode
# directories become barcode01 in tmp_files/barcodes/
prepare_input_layout () {
    if find "${DIR}" -mindepth 1 -maxdepth 1 -type d -name "barcode*" | grep -q .; then
        echo "Found barcode directories in ${DIR}"
        return 0
    fi
    if ! find "${DIR}" -maxdepth 1 -type f \( -name "*.fastq.gz" -o -name "*.fastq" \) | grep -q .; then
        echo "ERROR - No barcode directories or fastq files found in ${DIR} - Please fix your input data path." >&2
        echo "ABORTING" >&2
        exit 1
    fi
    echo "WARNING - Found fastq files in ${DIR} - Formatting tmp file structure. Please prefer a barcodeXX directory structure"

    local -r tmp_barcodes="tmp_files/barcodes"
    mkdir -p "${tmp_barcodes}/barcode01"
    local file
    for file in "${DIR}"/*.fastq.gz; do
        if [[ -e "${file}" ]]; then
            ln -sf "$(realpath "${file}")" "${tmp_barcodes}/barcode01"
        fi
    done

    # the snakefile reads *.fastq.gz only
    for file in "${DIR}"/*.fastq; do
        if [[ -e "${file}" ]]; then
            gzip -c "${file}" > "${tmp_barcodes}/barcode01/$(basename "${file}").gz"
        fi
    done

    if [[ ! -s "${METADATA}/metadata.csv" ]]; then
        echo "WARNING - No metadata file provided. NanoASV will generate a dummy one. Please prefer using your own metadata file"
        head -n2 "${MOCK_DIR}/metadata.csv" > "${tmp_barcodes}/metadata.csv"
    else
        cp "${METADATA}/metadata.csv" "${tmp_barcodes}/metadata.csv"
    fi

    DIR="${tmp_barcodes}"
    METADATA="${DIR}"
}

# usage: reject_metadata MESSAGE
reject_metadata () {
    echo "${1}" >&2
    cat "${NANOASV_PATH}/config/requirements.txt" >&2
    exit 1
}

# usage: check_metadata INPUT_DIR METADATA_DIR
# runs in a subshell (cd); called as "check_metadata ... || exit 1", so
# set -e is off: every check exits explicitly
check_metadata () (
    local input_path
    input_path="$(realpath "${1}")"
    cd "${2}" || exit 1

    [[ -s metadata.csv ]] || \
        { cowpy -e dead "Error : Please provide a metadata.csv" >&2 ; cat "${NANOASV_PATH}/config/requirements.txt" >&2 ; exit 1 ; }

    # at least 3 columns: row names and two data columns
    awk -F "," 'NR == 1 { exit NF > 2 ? 0 : 1}' metadata.csv || \
        reject_metadata 'ERROR: Check metadata.csv: it does not look like a csv file. Are you sure you are using coma to separate the fields? Do you have more than two columns?'

    awk -F "," 'NR == 1 { exit $1 == "" ? 0 : 1}' metadata.csv || \
        reject_metadata 'ERROR: First field of first line should be empty. Please check metadata.csv file structure.'

    awk 'END{ exit NR > 1 ? 0 : 1}' metadata.csv || \
        reject_metadata 'ERROR: metadata.csv: Missing header and/or data information. Too few lines.'

    local sample_name
    cut -f1 -d "," metadata.csv | \
        tail -n +2 | \
        while read -r sample_name ; do
            [[ -d "${input_path}/${sample_name}" ]] || \
                reject_metadata "ERROR, ${sample_name} sequences not found. Please check metadata.csv and barcodes directories"
        done || exit 1

    # the R step fails on a sample missing from metadata.csv
    local sample_dir
    for sample_dir in "${input_path}"/*/ ; do
        compgen -G "${sample_dir}*.fastq*" > /dev/null || continue
        sample_name="$(basename "${sample_dir}")"
        cut -f1 -d "," metadata.csv | \
            tail -n +2 | \
            grep -qxF "${sample_name}" || \
            reject_metadata "ERROR, ${sample_name} is not listed in metadata.csv. Please check metadata.csv and barcodes directories"
    done

    awk -F "," '{print NF}' metadata.csv | \
        sort -u | \
        awk 'END {exit NR == 1 ? 0 : 1}' || \
        reject_metadata 'ERROR: Check metadata.csv: not all the lines have the same number of columns'
)

# usage: parse_options "$@" (sets the option variables)
parse_options () {
    SNAKEMAKE_OPTIONS=()
    local key
    while [[ $# -gt 0 ]]; do
        key="$1"
        case "${key}" in
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
                echo 'NanoASV 1.2.3 - Conda-x-Snakemake - https://github.com/ImagoXV/NanoASV - Arthur Cousson, Frederic Mahe and Ulysse Guyet '
                exit
                ;;
            -h|--help)
                cat "${NANOASV_PATH}/config/help.txt"
                exit
                ;;
            --metadata)
                METADATA="$2"
                shift
                shift
                ;;
            --requirements)
                cat "${NANOASV_PATH}/config/requirements.txt"
                exit
                ;;
            --dry-run)
                DRY="n"
                shift
                ;;
            --mock)
                DIR="${MOCK_DIR}"
                OUT="${MOCK_OUT}"
                DATABASE="${MOCK_DATABASE}"
                NUM_PROCESSES="${MOCK_NUM_PROCESSES}"
                TMP_FILES=0
                shift
                ;;
            --model)
                MOD="$2"
                shift
                shift
                ;;
            --rerun-incomplete)
                SNAKEMAKE_OPTIONS+=(--rerun-incomplete)
                shift
                ;;
            --unlock)
                SNAKEMAKE_OPTIONS+=(--unlock)
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
                echo 'Unknown option:' "$1" >&2
                cat "${NANOASV_PATH}/config/help.txt" >&2
                exit 1
                ;;
        esac
    done

    # after the loop: options given before --dry-run are kept
    if [[ "${DRY}" == "n" ]]; then
        DIR="${DIR:-${MOCK_DIR}}"
        OUT="${OUT:-${MOCK_OUT}}"
        DATABASE="${DATABASE:-${MOCK_DATABASE}}"
    fi
}

# usage: set_defaults (fills in the options left out)
set_defaults () {
    QUAL="${QUAL:-$DEFAULT_QUAL}"
    MINL="${MINL:-$DEFAULT_MINL}"
    MAXL="${MAXL:-$DEFAULT_MAXL}"
    ID="${ID:-$DEFAULT_ID}"
    NUM_PROCESSES="${NUM_PROCESSES:-$DEFAULT_NUM_PROCESSES}"
    R_CLEANING="${R_CLEANING:-$DEFAULT_R_CLEANING}"
    SUBSAMPLING="${SUBSAMPLING:-$DEFAULT_SUBSAMPLING}"
    TREE="${TREE:-$DEFAULT_TREE}"
    METADATA="${METADATA:-${DIR}}"
    MINAB="${MINAB:-$DEFAULT_MINAB}"
    TMP_FILES="${TMP_FILES:-$DEFAULT_TMP_FILES}"
    MOD="${MOD:-$DEFAULT_MOD}"
    DATABASE="${DATABASE:-$DEFAULT_DATABASE}"
    FASTTREE_MOD="${FASTTREE_MOD:-$DEFAULT_FASTTREE_MOD}"
    SAMQ="${SAMQ:-$DEFAULT_SAMQ}"
}

check_subsampling () {
    if [[ ! "${SUBSAMPLING}" =~ ^[1-9][0-9]*$ ]]; then
        echo "ERROR: --subsampling must be a positive integer (number of reads per sample), not '${SUBSAMPLING}'." >&2
        exit 1
    fi
}

check_required_options () {
    if [[ -z "${DIR}" ]]; then
        cowpy -e dead "Error: -d needs an argument, I don't know where your sequences are." >&2
        cat "${NANOASV_PATH}/config/help.txt" >&2
        exit 1
    fi
    if [[ -z "${OUT}" ]]; then
        cowpy -e dead "Error: -o needs an argument. You don't want me to print to stdout" >&2
        exit 1
    fi
}

run_snakemake () {
    snakemake -"${DRY}"p -s "${NANOASV_PATH}/workflow/snakefile" "${SNAKEMAKE_OPTIONS[@]}" \
        --cores "${NUM_PROCESSES}" \
        --config \
            QUAL="${QUAL}" \
            MINL="${MINL}" \
            MAXL="${MAXL}" \
            ID="${ID}" \
            NUM_PROCESSES="${NUM_PROCESSES}" \
            R_CLEANING="${R_CLEANING}" \
            MINAB="${MINAB}" \
            SUBSAMPLING="$((SUBSAMPLING * LINES_PER_FASTQ_ENTRY))" \
            TREE="${TREE}" \
            INPUT_DIR="${DIR}" \
            OUT="${OUT}" \
            METADATA="${METADATA}" \
            DATABASE="${DATABASE}" \
            NANOASV_PATH="${NANOASV_PATH}" \
            MOD="${MOD}" \
            FASTTREE_MOD="${FASTTREE_MOD}" \
            SAMQ="${SAMQ}"
}

clean_up () {
    if [[ "${TMP_FILES}" -eq 0 ]]; then
        rm -r tmp_files
    fi

    if [[ "${DRY}" != "n" ]]; then
        tree "${OUT}"
    fi

    rm -rf tmp_files/tmp_seq_dir
}

main () {
    parse_options "$@"
    set_defaults
    check_subsampling
    mkdir -p tmp_files
    check_required_options

    gzip_fastq_in_barcode_dirs
    prepare_input_layout
    check_metadata "${DIR}" "${METADATA}" || exit 1

    check_minimap2_model
    check_fasttree_model
    check_reference_exists
    echo 'Checking reference format'
    check_reference_headers
    add_reference_ids
    linearize_reference

    run_snakemake
    clean_up
}

main "$@"
