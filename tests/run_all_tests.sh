#!/bin/bash
#
# Usage: bash tests/run_all_tests.sh (in an active NanoASV environment)

set -euo pipefail

TESTS_DIR="$(dirname "${BASH_SOURCE[0]}")"
readonly TESTS_DIR

for script in "${TESTS_DIR}"/{test_options,test_input}.sh ; do
    bash "${script}"
    echo
done
