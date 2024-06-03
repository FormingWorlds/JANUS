#!/bin/bash
# Usage: bash insert_star_spectrum.sh

read -p "Please enter Janus directory: " JANUS_DIR

# Set stellar spectrum to use
# This is the main parameter to change 
STAR_SPEC="$JANUS_DIR/src/janus/data/spectral_files/stellar_spectra/Sun_t4_4Ga_claire_12.txt"

# Set spectral file to use as basis
ORIG_FILE="$JANUS_DIR/src/janus/data/spectral_files/sp_b318_HITRAN_a16/sp_b318_HITRAN_a16_no_spectrum"

# New spectral file name
NEW_FILE="$JANUS_DIR/output/runtime_spectral_file"

# Insert stellar spectrum
cp ${ORIG_FILE} ${NEW_FILE}
cp ${ORIG_FILE}_k ${NEW_FILE}_k

prep_spec <<EOF
${NEW_FILE}
a
6
n
T
100 4000
100
2
n
${STAR_SPEC}
y
-1
EOF
