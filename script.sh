#!/bin/bash

# This script is the first docker version of NanoASV
# Authors : Arthur Cousson, Frederic Mahe
# 08/03/2023



# Manual entries - Arguments
# Set default values
DEFAULT_DIR="default/path/to/directory"
DEFAULT_OUT="default_output.txt"

# Read the arguments passed to the script
while [[ $# -gt 0 ]]; do
  key="$1"
  case $key in
    -d|--dir)
      DIR="$2"
      shift
      shift
      ;;
    -o|--out)
      OUT="$2"
      shift
      shift
      ;;
    *)
      echo "Unknown option: $1"
      shift
      ;;
  esac
done

# Assign default values if variables are empty
OUT="${OUT:-$DEFAULT_OUT}"

# Check if DIR is empty and no default value is provided
if [[ -z $DIR ]]; then
  echo "Error: -d needs an argument, I don't know where your sequences are."
  exit 1
fi



# Clean and create the Data directory
rm -vrf "$data_dir"
mkdir -vp "$data_dir"

# Create temporary directory

rm -vr ~/.tmp_NanoASV

date
echo Creating temporary directory at ~/.
mkdir ~/.tmp_NanoASV

# First step : 
