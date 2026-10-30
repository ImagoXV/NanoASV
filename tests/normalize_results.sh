#!/bin/bash
#
# Usage: normalize_results.sh RESULTS_DIR OUTPUT_DIR
#
# Copy the outputs of a NanoASV run into a form that can be compared
# with the golden files in tests/golden/: what changes from run to run
# (read order, tree branch lengths, unknown-cluster labels, binary
# .rdata) is normalized.

set -euo pipefail

readonly RESULTS="${1}"
readonly OUTPUT="${2}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
readonly SCRIPT_DIR
readonly UNKNOWN_CLUSTER_LABEL='[0-9a-f]{8}(-[0-9a-f]{4}){3}-[0-9a-f]{12}_Unknown_cluster_'
readonly NORMALIZED_UNKNOWN_CLUSTER_LABEL='UUID_Unknown_cluster_'
readonly NORMALIZE_LABELS="s/${UNKNOWN_CLUSTER_LABEL}/${NORMALIZED_UNKNOWN_CLUSTER_LABEL}/g"

mkdir -p "${OUTPUT}"/{ASV,Tax,CSV,Exact_affiliations,Phylogeny,Rdata}

(cd "${RESULTS}" && find . -type f | LC_ALL=C sort) > "${OUTPUT}/files.txt"

cp "${RESULTS}"/ASV/*.tsv "${OUTPUT}/ASV/"
cp "${RESULTS}"/Tax/*.csv "${OUTPUT}/Tax/"
for file in "${RESULTS}"/CSV/*.csv ; do
    sed -E "${NORMALIZE_LABELS}" \
        "${file}" > "${OUTPUT}/CSV/$(basename "${file}")"
done

if compgen -G "${RESULTS}/Unknown_clusters/*" > /dev/null ; then
    mkdir -p "${OUTPUT}/Unknown_clusters"
    for file in "${RESULTS}"/Unknown_clusters/* ; do
        sed -E "${NORMALIZE_LABELS}" \
            "${file}" > "${OUTPUT}/Unknown_clusters/$(basename "${file}")"
    done
fi

for file in "${RESULTS}"/Exact_affiliations/*.tsv ; do
    LC_ALL=C sort "${file}" > "${OUTPUT}/Exact_affiliations/$(basename "${file}")"
done

# keep the topology: drop branch lengths and support values
sed -E -e 's/:[-+.0-9eE]+//g ; s/\)[.0-9]+/)/g' \
       -e "${NORMALIZE_LABELS}" \
    "${RESULTS}/Phylogeny/ASV.tree" > "${OUTPUT}/Phylogeny/ASV.topology"

Rscript "${SCRIPT_DIR}/describe_rdata.r" \
    "${RESULTS}/Rdata/NanoASV.rdata" > "${OUTPUT}/Rdata/NanoASV.txt"
