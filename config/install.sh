#!/bin/bash
set -euo pipefail
INSTALL_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${INSTALL_DIR}/../"
VERSION="$(head -n1 environment.yml | cut -f2 -d " ")"
conda env create --file environment.yml -p "${VERSION}"
ACTIVATE_DIR="${INSTALL_DIR}/../${VERSION}/etc/conda/activate.d"
#ACTIVATE_DIR=$(conda env list | grep -w 'NanoASV' | awk '{print $2}')/etc/conda/activate.d
mkdir -p "${ACTIVATE_DIR}"
cp config/alias.sh "${ACTIVATE_DIR}/"
cp config/paths.sh "${ACTIVATE_DIR}/"
printf 'export NANOASV_PATH=%q\n' "$(pwd)" >> "${ACTIVATE_DIR}/paths.sh"
DEACTIVATE_DIR="${INSTALL_DIR}/../${VERSION}/etc/conda/deactivate.d"
#DEACTIVATE_DIR=$(conda env list | grep -w 'NanoASV' | awk '{print $2}')/etc/conda/deactivate.d
mkdir -p "${DEACTIVATE_DIR}"
cp config/unalias.sh "${DEACTIVATE_DIR}/"
chmod +x workflow/run.sh