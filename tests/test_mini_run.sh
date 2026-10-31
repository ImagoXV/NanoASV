#!/bin/bash
# shellcheck disable=SC2015
#
# Usage: bash tests/test_mini_run.sh (in an active NanoASV environment)
#
# Full runs on the dataset made by make_mini_dataset.py, where the fate
# of each read is known. Set NORMALIZED_MINI_RESULTS to a directory to
# keep the normalized results of the default run (golden files).

# shellcheck source=tests/common.sh
source "$(dirname "${BASH_SOURCE[0]}")/common.sh"
print_header "Mini run"

TESTS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
readonly TESTS_DIR
readonly DATASET="${TMP_ROOT}/mini"
readonly NORMALIZED="${NORMALIZED_MINI_RESULTS:-${TMP_ROOT}/normalized_mini}"

python3 -I "${TESTS_DIR}/make_mini_dataset.py" "${DATASET}"

# usage: mini_run INPUT_DIR OUTPUT_DIR [OPTION...]
mini_run () {
    local -r input_dir="${1}"
    local -r output_dir="${2}"
    shift 2
    "${NANOASV}" -d "${input_dir}" -o "${output_dir}" \
                 -db "${DATASET}/reference.fasta" "$@" > /dev/null 2>&1
}


## ------------------------------------------------------------ default run

DESCRIPTION="default run exits with status 0"
new_workdir
mini_run "${DATASET}" out && \
    success "${DESCRIPTION}" || \
        failure "${DESCRIPTION}"
RESULTS="${PWD}/out/Results"

DESCRIPTION="default run: normalized results match tests/golden/mini"
bash "${TESTS_DIR}/normalize_results.sh" "${RESULTS}" "${NORMALIZED}"
diff -r "${TESTS_DIR}/golden/mini" "${NORMALIZED}" && \
    success "${DESCRIPTION}" || \
        failure "${DESCRIPTION}"

DESCRIPTION="issue #154 --- reverse-strand reads are kept"
grep -qx "2 REV1\.1\.1500" "${RESULTS}/ASV/barcode01_abundance.tsv" && \
    success "${DESCRIPTION}" || \
        failure "${DESCRIPTION}"

# AB1.1.15 is a prefix of AB1.1.1500
for BARCODE in barcode01 barcode02 ; do
    DESCRIPTION="taxonomy of ${BARCODE}: exactly the references with reads"
    cmp -s <(cut -d " " -f 2 "${RESULTS}/ASV/${BARCODE}_abundance.tsv" | sort) \
           <(cut -d " " -f 1 "${RESULTS}/Tax/Taxonomy_${BARCODE}.csv" | sort) && \
        success "${DESCRIPTION}" || \
            failure "${DESCRIPTION}"
done

DESCRIPTION="5 identical unmatched reads form one unknown cluster"
tail -n +2 "${NORMALIZED}/Unknown_clusters/unknown_clusters.tsv" | \
    cmp -s - <(printf "UUID_Unknown_cluster_1\t5\n") && \
    success "${DESCRIPTION}" || \
        failure "${DESCRIPTION}"

DESCRIPTION="R cleaning removes Eukaryota, Mitochondria and Chloroplast"
grep -qE '^"(EUK|MITO|CHLO)\.' "${RESULTS}/CSV/Taxonomy-Abundance_table.csv" && \
    failure "${DESCRIPTION}" || \
        success "${DESCRIPTION}"

DESCRIPTION="R cleaning keeps taxa with missing ranks"
grep -q '^"SHORT\.1\.1500",' "${RESULTS}/CSV/Taxonomy-Abundance_table.csv" && \
    success "${DESCRIPTION}" || \
        failure "${DESCRIPTION}"


## ------------------------------------------------------- --no-r-cleaning

DESCRIPTION="--no-r-cleaning keeps Eukaryota, Mitochondria and Chloroplast"
new_workdir
mini_run "${DATASET}" out --no-r-cleaning
grep -cE '^"(EUK|MITO|CHLO)\.' out/Results/CSV/Taxonomy-Abundance_table.csv | \
    grep -qx "3" && \
    success "${DESCRIPTION}" || \
        failure "${DESCRIPTION}"


## ------------------------------------------------ exact ID matching

DESCRIPTION="the tree has no tip for a reference without reads (AB1.1.15 is a prefix of AB1.1.1500)"
new_workdir
cp -r "${DATASET}" input
rm -r input/barcode02
sed -i "/^barcode02,/d" input/metadata.csv
mini_run input out && \
    ! grep -qF "AB1.1.1500" out/Results/Phylogeny/ASV.tree && \
    success "${DESCRIPTION}" || \
        failure "${DESCRIPTION}"


## ------------------------------------------------ globs over tmp_files/

DESCRIPTION="tmp_files/ glob: a sample from an earlier run leaks into the results (known bug)"
new_workdir
cp -r "${DATASET}" with_barcode03
cp -r with_barcode03/barcode01 with_barcode03/barcode03
printf "barcode03,3,Mini\n" >> with_barcode03/metadata.csv
mini_run with_barcode03 out1
mini_run "${DATASET}" out2
[[ -e out2/Results/ASV/barcode03_abundance.tsv ]] && \
    success "${DESCRIPTION}" || \
        failure "${DESCRIPTION}"

DESCRIPTION="vsearch glob: unmatched reads of a sample not named barcode* are not clustered (known bug)"
new_workdir
cp -r "${DATASET}" input
mv input/barcode02 input/soilA
sed -i "s/^barcode02,/soilA,/" input/metadata.csv
mini_run input out && \
    ! compgen -G "out/Results/Unknown_clusters/*" > /dev/null && \
    success "${DESCRIPTION}" || \
        failure "${DESCRIPTION}"

exit_with_failure_count
