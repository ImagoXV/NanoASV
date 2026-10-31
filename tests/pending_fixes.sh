#!/bin/bash
# shellcheck disable=SC2015
#
# Usage: bash tests/pending_fixes.sh (in an active NanoASV environment)
#
# Known bugs, written as tests of the expected behaviour: they fail
# until the bug is fixed. Not run by run_all_tests.sh. When a bug is
# fixed, move its test to the matching test file.

# shellcheck source=tests/common.sh
source "$(dirname "${BASH_SOURCE[0]}")/common.sh"
print_header "Pending fixes"


DESCRIPTION="missing -d is reported when stdin is not a terminal (cowpy)"
new_workdir
"${NANOASV}" -o out 2>&1 | \
    grep -q -- "-d needs an argument" && \
    success "${DESCRIPTION}" || \
        failure "${DESCRIPTION}"

exit_with_failure_count
