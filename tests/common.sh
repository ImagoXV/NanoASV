# shellcheck shell=bash disable=SC2015
# Sourced by the test scripts. Requires an active NanoASV environment.

RED="\033[1;31m"
GREEN="\033[1;32m"
NO_COLOR="\033[0m"

failure () {
    printf "%bFAIL%b: %s\n" "${RED}" "${NO_COLOR}" "${1}"
    exit 1
}

success () {
    printf "%bPASS%b: %s\n" "${GREEN}" "${NO_COLOR}" "${1}"
}

print_header () {
    local line
    line="$(printf '%.0s-' {1..76})"
    printf "# %s %s\n" "${line:${#1}}" "${1}"
}

readonly NANOASV="${NANOASV_PATH}/workflow/run.sh"

TMP_ROOT="$(mktemp -d)"
readonly TMP_ROOT
trap 'rm -rf "${TMP_ROOT}"' EXIT

# run.sh calls conda init, which edits the shell startup files in ${HOME}
export HOME="${TMP_ROOT}/home"
mkdir -p "${HOME}"

# cowpy reads its message from stdin when stdin is not a terminal,
# and waits as long as stdin stays open
exec < /dev/null

# run.sh writes tmp_files/ in the current directory
new_workdir () {
    cd "$(mktemp -d -p "${TMP_ROOT}")" || exit 1
}

# usage: make_input_dir DIR
# two barcode directories with one read each, and a matching metadata.csv
make_input_dir () {
    local -r dir="${1}"
    local barcode
    for barcode in barcode01 barcode02 ; do
        mkdir -p "${dir}/${barcode}"
        printf "@r1\nACGT\n+\nIIII\n" | gzip > "${dir}/${barcode}/reads.fastq.gz"
    done
    printf ",ID,Origin\nbarcode01,1,A\nbarcode02,2,B\n" > "${dir}/metadata.csv"
}

# usage: dry_run_matches DESCRIPTION REGEX [OPTION...]
# run a dry run in a new working directory, search its output for REGEX
dry_run_matches () {
    local -r description="${1}"
    local -r regex="${2}"
    shift 2
    new_workdir
    "${NANOASV}" --dry-run "$@" 2>&1 | \
        grep -qE -- "${regex}" && \
        success "${description}" || \
            failure "${description}"
}
