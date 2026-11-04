#!/bin/bash
#
# Usage: read_counts.sh COUNTS_DIR UNKNOWN_CLUSTERS SAMPLE...
#
# Print the read-count table: one row per sample, one column per step.
# COUNTS_DIR/SAMPLE/STEP holds the number of reads after STEP.
# UNKNOWN_CLUSTERS (vsearch cluster table) is absent when no cluster is
# kept.

set -euo pipefail

readonly COUNTS_DIR="${1}"
readonly UNKNOWN_CLUSTERS="${2}"
shift 2
readonly -a STEPS=(input after_quality_filter after_adapter_trimming subsampled
                   unmapped low_mapq assigned)

# usage: reads_in_unknown_clusters SAMPLE
reads_in_unknown_clusters () {
    local -r sample="${1}"
    if [[ ! -s "${UNKNOWN_CLUSTERS}" ]] ; then
        echo 0
        return
    fi
    awk -F "\t" -v sample="${sample}" '
        NR == 1 {
            for (i = 2; i <= NF; i++) if ($i == sample) column = i
            next
        }
        column { sum += $column }
        END { print sum + 0 }' "${UNKNOWN_CLUSTERS}"
}

printf "sample"
printf "\t%s" "${STEPS[@]}" in_unknown_clusters
printf "\n"

for sample in "$@" ; do
    printf "%s" "${sample}"
    for step in "${STEPS[@]}" ; do
        printf "\t%s" "$(< "${COUNTS_DIR}/${sample}/${step}")"
    done
    printf "\t%s\n" "$(reads_in_unknown_clusters "${sample}")"
done
