#!/bin/bash
# shellcheck disable=SC2015
#
# Usage: bash tests/test_options.sh (in an active NanoASV environment)

# shellcheck source=tests/common.sh
source "$(dirname "${BASH_SOURCE[0]}")/common.sh"
print_header "Options"

DESCRIPTION="nanoasv is executable"
[[ -x "${NANOASV}" ]] && success "${DESCRIPTION}" || failure "${DESCRIPTION}"


## ------------------------------------------- help, version, requirements

for OPTION in -h --help ; do
    DESCRIPTION="${OPTION} prints config/help.txt"
    new_workdir
    "${NANOASV}" "${OPTION}" 2> /dev/null | \
        cmp -s - "${NANOASV_PATH}/config/help.txt" && \
        success "${DESCRIPTION}" || \
            failure "${DESCRIPTION}"

    DESCRIPTION="${OPTION} exits with status 0"
    new_workdir
    "${NANOASV}" "${OPTION}" > /dev/null 2>&1 && \
        success "${DESCRIPTION}" || \
            failure "${DESCRIPTION}"
done

for OPTION in -v --version ; do
    DESCRIPTION="${OPTION} prints the version number"
    new_workdir
    "${NANOASV}" "${OPTION}" 2> /dev/null | \
        grep -q "^NanoASV 1\.2\.2 " && \
        success "${DESCRIPTION}" || \
            failure "${DESCRIPTION}"
done

DESCRIPTION="nanoasv does not modify shell startup files"
new_workdir
"${NANOASV}" --version > /dev/null 2>&1
[[ ! -e "${HOME}/.bashrc" ]] && \
    success "${DESCRIPTION}" || \
        failure "${DESCRIPTION}"

DESCRIPTION="--requirements prints config/requirements.txt"
new_workdir
"${NANOASV}" --requirements 2> /dev/null | \
    cmp -s - "${NANOASV_PATH}/config/requirements.txt" && \
    success "${DESCRIPTION}" || \
        failure "${DESCRIPTION}"


## ------------------------------------------------ rejected command lines

DESCRIPTION="an unknown option is reported"
new_workdir
"${NANOASV}" --unknown 2> /dev/null | \
    grep -q "^Unknown option: --unknown$" && \
    success "${DESCRIPTION}" || \
        failure "${DESCRIPTION}"

DESCRIPTION="an unknown option exits with status 1"
new_workdir
"${NANOASV}" --unknown > /dev/null 2>&1 && \
    failure "${DESCRIPTION}" || \
        success "${DESCRIPTION}"

# help.txt documents --sam-quality, run.sh only accepts --sam-qual
DESCRIPTION="--sam-quality is an unknown option"
new_workdir
"${NANOASV}" --sam-quality 10 2> /dev/null | \
    grep -q "^Unknown option: --sam-quality$" && \
    success "${DESCRIPTION}" || \
        failure "${DESCRIPTION}"

DESCRIPTION="missing -d prints config/help.txt"
new_workdir
"${NANOASV}" -o out 2> /dev/null | \
    grep -q "Synopsis" && \
    success "${DESCRIPTION}" || \
        failure "${DESCRIPTION}"

DESCRIPTION="missing -d exits with status 1"
new_workdir
"${NANOASV}" -o out > /dev/null 2>&1 && \
    failure "${DESCRIPTION}" || \
        success "${DESCRIPTION}"

DESCRIPTION="missing -o exits with status 1"
new_workdir
make_input_dir input
"${NANOASV}" -d input > /dev/null 2>&1 && \
    failure "${DESCRIPTION}" || \
        success "${DESCRIPTION}"


## --------------------------------------------------------------- dry run

DESCRIPTION="--dry-run exits with status 0"
new_workdir
"${NANOASV}" --dry-run > /dev/null 2>&1 && \
    success "${DESCRIPTION}" || \
        failure "${DESCRIPTION}"

dry_run_matches "--dry-run plans all steps up to the R script" \
                "script\.r "


## --------------------------------------------------------- default values

dry_run_matches "default chopper parameters: -q 8 -l 1300 --maxlength 1700" \
                "chopper -q 8 -l 1300 --maxlength 1700 "

dry_run_matches "default subsampling: 50,000 reads (200,000 lines)" \
                "head -n 200000 "

dry_run_matches "default samtools mapping quality: 0" \
                "samtools view .* -q 0 "

dry_run_matches "default minimap2 model: map-ont" \
                "minimap2 -a -x map-ont "

dry_run_matches "default FastTree model: fastest" \
                "FastTree -nt -fastest "

dry_run_matches "default vsearch identity: 0.7" \
                "--id 0\.7 "

dry_run_matches "default minimum cluster abundance: 5" \
                "--minsize 5 "

dry_run_matches "default number of threads: 1" \
                "--threads 1 "

# mafft asks for 20 threads: Snakemake gives it at most --cores
dry_run_matches "default number of cores: 1 (mafft gets 1 thread)" \
                "mafft --thread 1 "

dry_run_matches "default R cleaning and tree: on (1 1)" \
                "script\.r [^ ]+ [^ ]+ 1 1 "


## ------------------------------------------------------ option wiring

for OPTION in -q --quality ; do
    dry_run_matches "${OPTION} is passed to chopper -q" \
                    "chopper -q 12 " "${OPTION}" 12
done

for OPTION in -l --minlength ; do
    dry_run_matches "${OPTION} is passed to chopper -l" \
                    "chopper -q 8 -l 1000 " "${OPTION}" 1000
done

for OPTION in -L --maxlength ; do
    dry_run_matches "${OPTION} is passed to chopper --maxlength" \
                    "--maxlength 2000 " "${OPTION}" 2000
done

dry_run_matches "--subsampling is converted to fastq lines (x 4)" \
                "head -n 40 " --subsampling 10

dry_run_matches "--sam-qual is passed to samtools view -q" \
                "samtools view .* -q 7 " --sam-qual 7

dry_run_matches "--model is passed to minimap2 -x" \
                "minimap2 -a -x asm5 " --model asm5

dry_run_matches "--fasttree-model is passed to FastTree" \
                "FastTree -nt -gtr " --fasttree-model gtr

for OPTION in -i --id-vsearch ; do
    dry_run_matches "${OPTION} is passed to vsearch --id" \
                    "--id 0\.9 " "${OPTION}" 0.9
done

for OPTION in -ab --minab ; do
    dry_run_matches "${OPTION} is passed to vsearch --minsize" \
                    "--minsize 7 " "${OPTION}" 7
done

for OPTION in -p --num-process ; do
    dry_run_matches "${OPTION} is passed to vsearch --threads" \
                    "--threads 3 " "${OPTION}" 3
    dry_run_matches "${OPTION} is passed to chopper --threads" \
                    "chopper .* --threads 3 " "${OPTION}" 3
    dry_run_matches "${OPTION} is passed to porechop --threads" \
                    "porechop .* --threads 3 " "${OPTION}" 3
    dry_run_matches "${OPTION} is passed to minimap2 -t (indexing)" \
                    "minimap2 -x map-ont -t 3 -d " "${OPTION}" 3
    dry_run_matches "${OPTION} is passed to minimap2 -t (mapping)" \
                    "minimap2 -a -x map-ont -t 3 " "${OPTION}" 3
    dry_run_matches "${OPTION} sets the number of cores (mafft gets 3 threads)" \
                    "mafft --thread 3 " "${OPTION}" 3
done

dry_run_matches "--no-r-cleaning is passed to the R script" \
                "script\.r [^ ]+ [^ ]+ 0 1 " --no-r-cleaning

dry_run_matches "--notree is passed to the R script" \
                "script\.r [^ ]+ [^ ]+ 1 0 " --notree

dry_run_matches "--notree skips the phylogeny" \
                "if \[ 0 -eq 1 \]" --notree

for OPTION in -d --dir ; do
    make_input_dir "${TMP_ROOT}/input_${OPTION}"
    dry_run_matches "${OPTION} is passed to the R script" \
                    "script\.r ${TMP_ROOT}/input_${OPTION} " \
                    "${OPTION}" "${TMP_ROOT}/input_${OPTION}"
done

for OPTION in -o --out ; do
    dry_run_matches "${OPTION} is passed to the R script" \
                    "script\.r [^ ]+ ${TMP_ROOT}/out " \
                    "${OPTION}" "${TMP_ROOT}/out"
done

for OPTION in -db --database ; do
    cp "${NANOASV_PATH}/config/MOCK/mock_references/complete-reference.fasta" \
       "${TMP_ROOT}/my_reference.fasta"
    dry_run_matches "${OPTION} is indexed by minimap2" \
                    "my_reference\.mmi ${TMP_ROOT}/my_reference\.fasta" \
                    "${OPTION}" "${TMP_ROOT}/my_reference.fasta"
done

mkdir "${TMP_ROOT}/metadata"
cp "${NANOASV_PATH}/config/MOCK/metadata.csv" "${TMP_ROOT}/metadata/"
dry_run_matches "--metadata is passed to the R script" \
                "script\.r .* ${TMP_ROOT}/metadata$" \
                --metadata "${TMP_ROOT}/metadata"

exit 0
