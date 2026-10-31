# shellcheck shell=bash disable=SC2015
# Sourced by the test scripts. Requires an active NanoASV environment.

RED="\033[1;31m"
GREEN="\033[1;32m"
NO_COLOR="\033[0m"

FAILURES=0

failure () {
    printf "%bFAIL%b: %s\n" "${RED}" "${NO_COLOR}" "${1}"
    FAILURES=$((FAILURES + 1))
}

success () {
    printf "%bPASS%b: %s\n" "${GREEN}" "${NO_COLOR}" "${1}"
}

# usage: exit_with_failure_count (last line of a test script)
exit_with_failure_count () {
    if (( FAILURES > 0 )) ; then
        printf "%b%d test(s) failed%b\n" "${RED}" "${FAILURES}" "${NO_COLOR}"
        exit 1
    fi
    exit 0
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

# usage: exits_with STATUS COMMAND [ARG...]
# run COMMAND, output discarded; succeed if its exit status is STATUS
exits_with () {
    local -r expected_status="${1}"
    shift
    local status=0
    "$@" > /dev/null 2>&1 || status=$?
    (( status == expected_status ))
}

# usage: rule_jobs RULE < DRY_RUN_OUTPUT
# keep the jobs of RULE: each job starts with a timestamp line "[...]"
rule_jobs () {
    awk -v rule_line="rule ${1}:" '
        /^\[/ { in_rule = 0 }
        $0 == rule_line { in_rule = 1 }
        in_rule'
}

# usage: dry_run_matches DESCRIPTION RULE REGEX [OPTION...]
# run a dry run in a new working directory, search the jobs of RULE for REGEX
dry_run_matches () {
    local -r description="${1}"
    local -r rule="${2}"
    local -r regex="${3}"
    shift 3
    new_workdir
    "${NANOASV}" --dry-run "$@" 2>&1 | \
        rule_jobs "${rule}" | \
        grep -qE -- "${regex}" && \
        success "${description}" || \
            failure "${description}"
}
