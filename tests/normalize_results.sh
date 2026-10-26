#!/bin/bash
#
# Usage: normalize_results.sh RESULTS_DIR OUTPUT_DIR
#
# Copy the outputs of a NanoASV run into a form that can be compared
# with the golden files in tests/golden/: what changes from run to run
# (read order, tree branch lengths, binary .rdata) is normalized.

set -euo pipefail

readonly RESULTS="${1}"
readonly OUTPUT="${2}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
readonly SCRIPT_DIR

mkdir -p "${OUTPUT}"/{ASV,Tax,CSV,Exact_affiliations,Phylogeny,Rdata}

(cd "${RESULTS}" && find . -type f | LC_ALL=C sort) > "${OUTPUT}/files.txt"

cp "${RESULTS}"/ASV/*.tsv "${OUTPUT}/ASV/"
cp "${RESULTS}"/Tax/*.csv "${OUTPUT}/Tax/"
cp "${RESULTS}"/CSV/*.csv "${OUTPUT}/CSV/"

for file in "${RESULTS}"/Exact_affiliations/*.tsv ; do
    LC_ALL=C sort "${file}" > "${OUTPUT}/Exact_affiliations/$(basename "${file}")"
done

# keep the topology: drop branch lengths and support values
sed -E 's/:[-+.0-9eE]+//g ; s/\)[.0-9]+/)/g' \
    "${RESULTS}/Phylogeny/ASV.tree" > "${OUTPUT}/Phylogeny/ASV.topology"

Rscript "${SCRIPT_DIR}/describe_rdata.r" \
    "${RESULTS}/Rdata/NanoASV.rdata" > "${OUTPUT}/Rdata/NanoASV.txt"
