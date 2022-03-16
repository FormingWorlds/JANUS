# Usage: bash insert_star_spectrum.sh

# Define path and stellar spectrum to use
STAR_PATH=/Users/timlichtenberg/git/proteus/atm_rad_conv/spectral_files/stellar_spectra/

# Options: Sun_t456myr_kurucz_95, F2V_hd128167, M45_ADLeo
STAR_NAME=Sun_t456myr_kurucz_95

# File name ending
STAR_FILE=.txt

# New spectral file name
NEW_FILE=sp_b318_HITRAN_a16_${STAR_NAME}

cp sp_b318_HITRAN_a16_no_spectrum ${NEW_FILE}
cp sp_b318_HITRAN_a16_no_spectrum_k ${NEW_FILE}_k

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
${STAR_PATH}${STAR_NAME}${STAR_FILE}
y
-1
EOF
