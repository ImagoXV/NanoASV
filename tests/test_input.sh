#!/bin/bash
# shellcheck disable=SC2015
#
# Usage: bash tests/test_input.sh (in an active NanoASV environment)

# shellcheck source=tests/common.sh
source "$(dirname "${BASH_SOURCE[0]}")/common.sh"
print_header "Input checks"

readonly MOCK_REFERENCE="${NANOASV_PATH}/config/MOCK/mock_references/complete-reference.fasta"


## ------------------------------------------------------- input directory

DESCRIPTION="a missing input directory is reported"
new_workdir
"${NANOASV}" -d missing -o out 2> /dev/null | \
    grep -q "^ERROR - No barcode directories or fastq files found in missing" && \
    success "${DESCRIPTION}" || \
        failure "${DESCRIPTION}"

DESCRIPTION="a missing input directory exits with status 1"
new_workdir
"${NANOASV}" -d missing -o out > /dev/null 2>&1 && \
    failure "${DESCRIPTION}" || \
        success "${DESCRIPTION}"

DESCRIPTION="an empty input directory is reported"
new_workdir
mkdir empty
"${NANOASV}" -d empty -o out 2> /dev/null | \
    grep -q "^ERROR - No barcode directories or fastq files found in empty" && \
    success "${DESCRIPTION}" || \
        failure "${DESCRIPTION}"

DESCRIPTION="barcode directories are found"
new_workdir
make_input_dir input
"${NANOASV}" --dry-run -d input 2> /dev/null | \
    grep -q "^Found barcode directories in input$" && \
    success "${DESCRIPTION}" || \
        failure "${DESCRIPTION}"

DESCRIPTION="fastq files without barcode directories are moved to barcode01"
new_workdir
mkdir input
printf "@r1\nACGT\n+\nIIII\n" | gzip > input/reads.fastq.gz
"${NANOASV}" --dry-run -d input > /dev/null 2>&1
[[ -L tmp_files/barcodes/barcode01/reads.fastq.gz ]] && \
    success "${DESCRIPTION}" || \
        failure "${DESCRIPTION}"

DESCRIPTION="fastq files without barcode directories trigger a warning"
new_workdir
mkdir input
printf "@r1\nACGT\n+\nIIII\n" | gzip > input/reads.fastq.gz
"${NANOASV}" --dry-run -d input 2> /dev/null | \
    grep -q "^WARNING - Found fastq files in input" && \
    success "${DESCRIPTION}" || \
        failure "${DESCRIPTION}"

DESCRIPTION="without metadata.csv, a dummy one is made from the mock metadata"
new_workdir
mkdir input
printf "@r1\nACGT\n+\nIIII\n" | gzip > input/reads.fastq.gz
"${NANOASV}" --dry-run -d input > /dev/null 2>&1
head -n 2 "${NANOASV_PATH}/config/MOCK/metadata.csv" | \
    cmp -s - tmp_files/barcodes/metadata.csv && \
    success "${DESCRIPTION}" || \
        failure "${DESCRIPTION}"

DESCRIPTION="uncompressed fastq files are reported"
new_workdir
make_input_dir input
gunzip input/barcode01/reads.fastq.gz
"${NANOASV}" --dry-run -d input 2> /dev/null | \
    grep -q "^NanoASV found flat fastq files. Compressing" && \
    success "${DESCRIPTION}" || \
        failure "${DESCRIPTION}"

# tmp_files/tmp_seq_dir/ is removed when run.sh ends
DESCRIPTION="uncompressed fastq files are compressed and used as input"
new_workdir
make_input_dir input
gunzip input/barcode01/reads.fastq.gz
"${NANOASV}" --dry-run -d input 2> /dev/null | \
    grep -q "^Found barcode directories in tmp_files/tmp_seq_dir$" && \
    success "${DESCRIPTION}" || \
        failure "${DESCRIPTION}"


## ----------------------------------------------------------- metadata.csv

# usage: metadata_rejected DESCRIPTION MESSAGE_REGEX METADATA_CONTENT
metadata_rejected () {
    local -r description="${1}"
    local -r regex="${2}"
    local -r content="${3}"
    new_workdir
    make_input_dir input
    printf "%b" "${content}" > input/metadata.csv
    "${NANOASV}" --dry-run -d input 2> /dev/null | \
        grep -qE -- "${regex}" && \
        success "${description}" || \
            failure "${description}"

    new_workdir
    make_input_dir input
    printf "%b" "${content}" > input/metadata.csv
    "${NANOASV}" --dry-run -d input > /dev/null 2>&1 && \
        failure "${description} (exit status 1)" || \
            success "${description} (exit status 1)"
}

DESCRIPTION="a missing metadata.csv is rejected"
new_workdir
make_input_dir input
rm input/metadata.csv
"${NANOASV}" --dry-run -d input > /dev/null 2>&1 && \
    failure "${DESCRIPTION}" || \
        success "${DESCRIPTION}"

DESCRIPTION="a missing metadata.csv prints config/requirements.txt"
new_workdir
make_input_dir input
rm input/metadata.csv
"${NANOASV}" --dry-run -d input 2> /dev/null | \
    grep -qF "$(head -n 1 "${NANOASV_PATH}/config/requirements.txt")" && \
    success "${DESCRIPTION}" || \
        failure "${DESCRIPTION}"

metadata_rejected "metadata.csv with less than 3 columns is rejected" \
                  "does not look like a csv file" \
                  ",ID\nbarcode01,1\nbarcode02,2\n"

metadata_rejected "metadata.csv with tabs instead of commas is rejected" \
                  "does not look like a csv file" \
                  "\tID\tOrigin\nbarcode01\t1\tA\nbarcode02\t2\tB\n"

metadata_rejected "metadata.csv with a non-empty first header field is rejected" \
                  "First field of first line should be empty" \
                  "sample,ID,Origin\nbarcode01,1,A\nbarcode02,2,B\n"

metadata_rejected "metadata.csv with a header only is rejected" \
                  "Too few lines" \
                  ",ID,Origin\n"

metadata_rejected "metadata.csv listing a missing barcode is rejected" \
                  "^ERROR, barcode03 sequences not found" \
                  ",ID,Origin\nbarcode01,1,A\nbarcode03,3,C\n"

metadata_rejected "metadata.csv with uneven numbers of columns is rejected" \
                  "not all the lines have the same number of columns" \
                  ",ID,Origin\nbarcode01,1,A\nbarcode02,2,B,extra\n"

DESCRIPTION="metadata.csv may list a subset of the barcodes"
new_workdir
make_input_dir input
printf ",ID,Origin\nbarcode01,1,A\n" > input/metadata.csv
"${NANOASV}" --dry-run -d input > /dev/null 2>&1 && \
    success "${DESCRIPTION}" || \
        failure "${DESCRIPTION}"


## --------------------------------------------------------- minimap2 model

DESCRIPTION="an unknown minimap2 model is reported"
new_workdir
"${NANOASV}" --dry-run --model unknown 2> /dev/null | \
    grep -q "^ERROR: --model invalid specification" && \
    success "${DESCRIPTION}" || \
        failure "${DESCRIPTION}"

DESCRIPTION="an unknown minimap2 model exits with status 1"
new_workdir
"${NANOASV}" --dry-run --model unknown > /dev/null 2>&1 && \
    failure "${DESCRIPTION}" || \
        success "${DESCRIPTION}"


## ------------------------------------------------------------- reference

DESCRIPTION="a missing reference file is reported"
new_workdir
"${NANOASV}" --dry-run -db missing.fasta 2> /dev/null | \
    grep -q "^No reference file found at missing.fasta" && \
    success "${DESCRIPTION}" || \
        failure "${DESCRIPTION}"

DESCRIPTION="a missing reference file exits with status 1"
new_workdir
"${NANOASV}" --dry-run -db missing.fasta > /dev/null 2>&1 && \
    failure "${DESCRIPTION}" || \
        success "${DESCRIPTION}"

DESCRIPTION="an empty reference file is rejected"
new_workdir
touch empty.fasta
"${NANOASV}" --dry-run -db empty.fasta > /dev/null 2>&1 && \
    failure "${DESCRIPTION}" || \
        success "${DESCRIPTION}"

DESCRIPTION="reference headers without IDs are reported"
new_workdir
printf ">Bacteria;Phylum;Class;Order;Family;Genus;Species\nACGT\n" > ref.fasta
"${NANOASV}" --dry-run -db ref.fasta 2> /dev/null | \
    grep -q "^Invalid reference file fasta header" && \
    success "${DESCRIPTION}" || \
        failure "${DESCRIPTION}"

DESCRIPTION="reference headers without IDs get a uuid"
new_workdir
printf ">Bacteria;Phylum;Class;Order;Family;Genus;Species\nACGT\n" > ref.fasta
"${NANOASV}" --dry-run -db ref.fasta > /dev/null 2>&1
grep -qE "^>[0-9a-f]{8}(-[0-9a-f]{4}){3}-[0-9a-f]{12} Bacteria;Phylum;" \
     tmp_files/CORRECT_HEADER_reference.fasta && \
    success "${DESCRIPTION}" || \
        failure "${DESCRIPTION}"

DESCRIPTION="reference headers without IDs get unique uuids"
new_workdir
awk 'BEGIN { for (i = 1; i <= 1000; i++) print ">Bacteria;P;C;O;F;G;S\nACGT" }' > ref.fasta
"${NANOASV}" --dry-run -db ref.fasta > /dev/null 2>&1
grep -o "^>[^ ]*" tmp_files/CORRECT_HEADER_reference.fasta | \
    sort -u | \
    wc -l | \
    grep -qx "1000" && \
    success "${DESCRIPTION}" || \
        failure "${DESCRIPTION}"

DESCRIPTION="reference headers without IDs: a file from an earlier run is replaced"
new_workdir
printf ">Bacteria;Phylum;Class;Order;Family;Genus;Species\nACGT\n" > ref.fasta
"${NANOASV}" --dry-run -db ref.fasta > /dev/null 2>&1
"${NANOASV}" --dry-run -db ref.fasta > /dev/null 2>&1
grep -c "^>" tmp_files/CORRECT_HEADER_reference.fasta | \
    grep -qx "1" && \
    success "${DESCRIPTION}" || \
        failure "${DESCRIPTION}"

DESCRIPTION="headers of a gzipped reference without IDs get a uuid"
new_workdir
printf ">Bacteria;Phylum;Class;Order;Family;Genus;Species\nACGT\n" | gzip > ref.fasta.gz
"${NANOASV}" --dry-run -db ref.fasta.gz > /dev/null 2>&1
grep -qE "^>[0-9a-f]{8}(-[0-9a-f]{4}){3}-[0-9a-f]{12} Bacteria;Phylum;" \
     tmp_files/CORRECT_HEADER_reference.fasta && \
    success "${DESCRIPTION}" || \
        failure "${DESCRIPTION}"

DESCRIPTION="headers of a gzipped reference are read after decompression"
new_workdir
printf ">s1 Bacteria;P;C;O;F;G;S\nACGT\n" | gzip > ref.fasta.gz
"${NANOASV}" --dry-run -db ref.fasta.gz > /dev/null 2>&1 && \
    success "${DESCRIPTION}" || \
        failure "${DESCRIPTION}"

DESCRIPTION="reference headers with IDs are kept"
new_workdir
"${NANOASV}" --dry-run -db "${MOCK_REFERENCE}" > /dev/null 2>&1
[[ ! -e tmp_files/CORRECT_HEADER_reference.fasta ]] && \
    success "${DESCRIPTION}" || \
        failure "${DESCRIPTION}"

DESCRIPTION="an interleaved reference is reported"
new_workdir
printf ">s1 Bacteria;P;C;O;F;G;S\nACGT\nACGT\n" > ref.fasta
"${NANOASV}" --dry-run -db ref.fasta 2> /dev/null | \
    grep -q "^Interleaved reference fasta file - " && \
    success "${DESCRIPTION}" || \
        failure "${DESCRIPTION}"

DESCRIPTION="an uncompressed interleaved reference is linearized"
new_workdir
printf ">s1 Bacteria;P;C;O;F;G;S\nACGT\nACGT\n" > ref.fasta
"${NANOASV}" --dry-run -db ref.fasta > /dev/null 2>&1
printf ">s1 Bacteria;P;C;O;F;G;S\nACGTACGT\n" | \
    cmp -s - tmp_files/SINGLELINE_reference.fasta && \
    success "${DESCRIPTION}" || \
        failure "${DESCRIPTION}"

DESCRIPTION="a gzipped interleaved reference is linearized"
new_workdir
printf ">s1 Bacteria;P;C;O;F;G;S\nACGT\nACGT\n" | gzip > ref.fasta.gz
"${NANOASV}" --dry-run -db ref.fasta.gz > /dev/null 2>&1
printf ">s1 Bacteria;P;C;O;F;G;S\nACGTACGT\n" | \
    cmp -s - tmp_files/SINGLELINE_reference.fasta && \
    success "${DESCRIPTION}" || \
        failure "${DESCRIPTION}"


exit 0
