#!/bin/bash
INSTALL_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd $INSTALL_DIR
mkdir conda
conda env create --file environment.yml -p conda
ACTIVATE_DIR=$(conda env list | grep -w 'NanoASV' | awk '{print $2}')/etc/conda/activate.d
cp config/alias.sh $ACTIVATE_DIR/
cp config/paths.sh $ACTIVATE_DIR/
echo "export NANOASV_PATH=$(pwd)" >> $ACTIVATE_DIR/paths.sh
DEACTIVATE_DIR=$(conda env list | grep -w 'NanoASV' | awk '{print $2}')/etc/conda/deactivate.d
cp config/unalias.sh $DEACTIVATE_DIR/
chmod +x workflow/run.sh