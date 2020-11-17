echo '###############   NH3 lines – HITRAN    ###############'
Ccorr_k -F pt_grid_lbl.dat -D /network/group/aopp/planetary/RTP015_LICHTENBERG_20GEDGE/spectral_files/dat_hitran/nh3_hitran.par -R 1 318 -c 3000.0 -i 0.1 -l 11 1.0e1 -t 1.0e-2 -k -s sp_b318_HITRAN_a16 +p -lk -o nh3_o -m nh3_m -L nh3_lbl.nc -np 16
prep_spec <<EOF
sp_b318_HITRAN_a16
a
5
y
nh3_o
-1
EOF
touch temp/done_nh3_lbl
rsync -av sp_b318_HITRAN_a16 temp/sp_b318_HITRAN_a16_11
rsync -av sp_b318_HITRAN_a16_k temp/sp_b318_HITRAN_a16_11_k
echo '###############   HNO3 lines – HITRAN    ###############'
Ccorr_k -F pt_grid_lbl.dat -D /network/group/aopp/planetary/RTP015_LICHTENBERG_20GEDGE/spectral_files/dat_hitran/hno3_hitran.par -R 1 318 -c 3000.0 -i 0.1 -l 12 1.0e1 -t 1.0e-2 -k -s sp_b318_HITRAN_a16 +p -lk -o hno3_o -m hno3_m -L hno3_lbl.nc -np 16
prep_spec <<EOF
sp_b318_HITRAN_a16
a
5
y
hno3_o
-1
EOF
touch temp/done_hno3_lbl
rsync -av sp_b318_HITRAN_a16 temp/sp_b318_HITRAN_a16_12
rsync -av sp_b318_HITRAN_a16_k temp/sp_b318_HITRAN_a16_12_k
echo '###############   N2 lines – HITRAN    ###############'
Ccorr_k -F pt_grid_lbl.dat -D /network/group/aopp/planetary/RTP015_LICHTENBERG_20GEDGE/spectral_files/dat_hitran/n2_hitran.par -R 1 318 -c 3000.0 -i 0.1 -l 13 1.0e1 -t 1.0e-2 -k -s sp_b318_HITRAN_a16 +p -lk -o n2_o -m n2_m -L n2_lbl.nc -np 16
prep_spec <<EOF
sp_b318_HITRAN_a16
a
5
y
n2_o
-1
EOF
touch temp/done_n2_lbl
rsync -av sp_b318_HITRAN_a16 temp/sp_b318_HITRAN_a16_13
rsync -av sp_b318_HITRAN_a16_k temp/sp_b318_HITRAN_a16_13_k
echo '###############   OCS lines – HITRAN    ###############'
Ccorr_k -F pt_grid_lbl.dat -D /network/group/aopp/planetary/RTP015_LICHTENBERG_20GEDGE/spectral_files/dat_hitran/ocs_hitran.par -R 1 318 -c 3000.0 -i 0.1 -l 25 1.0e1 -t 1.0e-2 -k -s sp_b318_HITRAN_a16 +p -lk -o ocs_o -m ocs_m -L ocs_lbl.nc -np 16
prep_spec <<EOF
sp_b318_HITRAN_a16
a
5
y
ocs_o
-1
EOF
touch temp/done_ocs_lbl
rsync -av sp_b318_HITRAN_a16 temp/sp_b318_HITRAN_a16_14
rsync -av sp_b318_HITRAN_a16_k temp/sp_b318_HITRAN_a16_14_k
echo '###############   N2-N2 CIA   ###############'
Ccorr_k -F pt_grid_cia.dat -CIA /network/group/aopp/planetary/RTP015_LICHTENBERG_20GEDGE/spectral_files/dat_continua/N2-N2_2018_fixed.cia -R 1 318 -i 0.1 -ct 13 13 1.0e2 -t 1.0e-2 -s sp_b318_HITRAN_a16 +p -lk -o n2n2_o -m n2n2_m -L n2n2_lbl.nc -np 16
prep_spec <<EOF
sp_b318_HITRAN_a16
a
19
y
n2n2_o
-1
EOF
sed -i -e 's/\*\*\*\*\*/0/g' sp_b318_HITRAN_a16
sed -i -e 's/\*\*\*//g' sp_b318_HITRAN_a16
sed -i -e 's/            NaN/0.000000000E+00/g' sp_b318_HITRAN_a16
touch temp/done_n2n2_cia
rsync -av sp_b318_HITRAN_a16 temp/sp_b318_HITRAN_a16_15
rsync -av sp_b318_HITRAN_a16_k temp/sp_b318_HITRAN_a16_15_k
echo '###############   O2-O2 CIA   ###############'
Ccorr_k -F pt_grid_cia.dat -CIA /network/group/aopp/planetary/RTP015_LICHTENBERG_20GEDGE/spectral_files/dat_continua/O2-O2_2018b_fixed.cia -R 1 318 -i 0.1 -ct 7 7 1.0e2 -t 1.0e-2 -s sp_b318_HITRAN_a16 +p -lk -o o2o2_o -m o2o2_m -L o2o2_lbl.nc -np 16
prep_spec <<EOF
sp_b318_HITRAN_a16
a
19
y
o2o2_o
-1
EOF
sed -i -e 's/\*\*\*\*\*/0/g' sp_b318_HITRAN_a16
sed -i -e 's/\*\*\*//g' sp_b318_HITRAN_a16
sed -i -e 's/            NaN/0.000000000E+00/g' sp_b318_HITRAN_a16
touch temp/done_o2o2_cia
rsync -av sp_b318_HITRAN_a16 temp/sp_b318_HITRAN_a16_16
rsync -av sp_b318_HITRAN_a16_k temp/sp_b318_HITRAN_a16_16_k
echo '###############   H2-H2 CIA   ###############'
Ccorr_k -F pt_grid_cia.dat -CIA /network/group/aopp/planetary/RTP015_LICHTENBERG_20GEDGE/spectral_files/dat_continua/H2-H2_CIA_Borysow_combined.cia -R 1 318 -i 0.1 -ct 23 23 1.0e2 -t 1.0e-2 -s sp_b318_HITRAN_a16 +p -lk -o h2h2_o -m h2h2_m -L h2h2_lbl.nc -np 16
prep_spec <<EOF
sp_b318_HITRAN_a16
a
19
y
h2h2_o
-1
EOF
sed -i -e 's/\*\*\*\*\*/0/g' sp_b318_HITRAN_a16
sed -i -e 's/\*\*\*//g' sp_b318_HITRAN_a16
sed -i -e 's/            NaN/0.000000000E+00/g' sp_b318_HITRAN_a16
touch temp/done_h2h2_cia
rsync -av sp_b318_HITRAN_a16 temp/sp_b318_HITRAN_a16_17
rsync -av sp_b318_HITRAN_a16_k temp/sp_b318_HITRAN_a16_17_k
echo '###############   CO2-CO2 CIA   ###############'
Ccorr_k -F pt_grid_cia.dat -CIA /network/group/aopp/planetary/RTP015_LICHTENBERG_20GEDGE/spectral_files/dat_continua/CO2-CO2_2018_fixed.cia -R 1 318 -i 0.1 -ct 2 2 1.0e2 -t 1.0e-2 -s sp_b318_HITRAN_a16 +p -lk -o co2co2_o -m co2co2_m -L co2co2_lbl.nc -np 16
prep_spec <<EOF
sp_b318_HITRAN_a16
a
19
y
co2co2_o
-1
EOF
sed -i -e 's/\*\*\*\*\*/0/g' sp_b318_HITRAN_a16
sed -i -e 's/\*\*\*//g' sp_b318_HITRAN_a16
sed -i -e 's/            NaN/0.000000000E+00/g' sp_b318_HITRAN_a16
touch temp/done_co2co2_cia
rsync -av sp_b318_HITRAN_a16 temp/sp_b318_HITRAN_a16_18
rsync -av sp_b318_HITRAN_a16_k temp/sp_b318_HITRAN_a16_18_k
echo '###############   CH4-CH4 CIA   ###############'
Ccorr_k -F pt_grid_cia.dat -CIA /network/group/aopp/planetary/RTP015_LICHTENBERG_20GEDGE/spectral_files/dat_continua/CH4-CH4_2011.cia -R 1 318 -i 0.1 -ct 6 6 1.0e2 -t 1.0e-2 -s sp_b318_HITRAN_a16 +p -lk -o ch4ch4_o -m ch4ch4_m -L ch4ch4_lbl.nc -np 16
prep_spec <<EOF
sp_b318_HITRAN_a16
a
19
y
ch4ch4_o
-1
EOF
sed -i -e 's/\*\*\*\*\*/0/g' sp_b318_HITRAN_a16
sed -i -e 's/\*\*\*//g' sp_b318_HITRAN_a16
sed -i -e 's/            NaN/0.000000000E+00/g' sp_b318_HITRAN_a16
touch temp/done_ch4ch4_cia
rsync -av sp_b318_HITRAN_a16 temp/sp_b318_HITRAN_a16_19
rsync -av sp_b318_HITRAN_a16_k temp/sp_b318_HITRAN_a16_19_k
echo '###############   H2O-H2O MT_CKD   ###############'
Ccorr_k -F pt_grid_cia.dat -R 1 318 -c 2500.0 -i 0.1 -ct 1 1 10.0 -t 1.0e-3 -e /network/group/aopp/planetary/RTP015_LICHTENBERG_20GEDGE/spectral_files/dat_continua/mt_ckd_v3.0_s296 /network/group/aopp/planetary/RTP015_LICHTENBERG_20GEDGE/spectral_files/dat_continua/mt_ckd_v3.0_s260 -k -s sp_b318_HITRAN_a16 +p -lk -o h2o-h2o_l318_1-318c -m h2o-h2o_l318_1-318cm -L h2o-h2o_lbl_lw.nc -lw h2o_l318_1-318map.nc -np 16
prep_spec <<EOF
sp_b318_HITRAN_a16
a
19
y
h2o-h2o_l318_1-318c
-1
EOF
sed -i -e 's/\*\*\*\*\*/0/g' sp_b318_HITRAN_a16
sed -i -e 's/\*\*\*//g' sp_b318_HITRAN_a16
sed -i -e 's/            NaN/0.000000000E+00/g' sp_b318_HITRAN_a16
touch temp/done_h2o_mtckd
rsync -av sp_b318_HITRAN_a16 temp/sp_b318_HITRAN_a16_20
rsync -av sp_b318_HITRAN_a16_k temp/sp_b318_HITRAN_a16_20_k
echo '###############   CO2-CH4 CIA   ###############'
Ccorr_k -F pt_grid_cia.dat -CIA /network/group/aopp/planetary/RTP015_LICHTENBERG_20GEDGE/spectral_files/dat_continua/CO2-CH4_2018.cia -R 1 318 -i 0.1 -ct 2 6 1.0e2 -t 1.0e-2 -s sp_b318_HITRAN_a16 +p -lk -o co2ch4_o -m co2ch4_m -L co2ch4_lbl.nc -np 16
prep_spec <<EOF
sp_b318_HITRAN_a16
a
19
y
co2ch4_o
-1
EOF
sed -i -e 's/\*\*\*\*\*/0/g' sp_b318_HITRAN_a16
sed -i -e 's/\*\*\*//g' sp_b318_HITRAN_a16
sed -i -e 's/            NaN/0.000000000E+00/g' sp_b318_HITRAN_a16
sed -i -e 's/          51/           0/g' sp_b318_HITRAN_a16
touch temp/done_co2ch4_cia
rsync -av sp_b318_HITRAN_a16 temp/sp_b318_HITRAN_a16_21
rsync -av sp_b318_HITRAN_a16_k temp/sp_b318_HITRAN_a16_21_k
echo '###############   CO2-H2 CIA   ###############'
Ccorr_k -F pt_grid_cia.dat -CIA /network/group/aopp/planetary/RTP015_LICHTENBERG_20GEDGE/spectral_files/dat_continua/CO2-H2_2018.cia -R 1 318 -i 0.1 -ct 2 23 1.0e2 -t 1.0e-2 -s sp_b318_HITRAN_a16 +p -lk -o co2h2_o -m co2h2_m -L co2h2_lbl.nc -np 16
prep_spec <<EOF
sp_b318_HITRAN_a16
a
19
y
co2h2_o
-1
EOF
sed -i -e 's/\*\*\*\*\*/0/g' sp_b318_HITRAN_a16
sed -i -e 's/\*\*\*//g' sp_b318_HITRAN_a16
sed -i -e 's/            NaN/0.000000000E+00/g' sp_b318_HITRAN_a16
sed -i -e 's/          51/           0/g' sp_b318_HITRAN_a16
touch temp/done_co2h2_cia
rsync -av sp_b318_HITRAN_a16 temp/sp_b318_HITRAN_a16_22
rsync -av sp_b318_HITRAN_a16_k temp/sp_b318_HITRAN_a16_22_k
echo '###############   CO2-He CIA   ###############'
Ccorr_k -F pt_grid_cia.dat -CIA /network/group/aopp/planetary/RTP015_LICHTENBERG_20GEDGE/spectral_files/dat_continua/CO2-He_2018.cia -R 1 318 -i 0.1 -ct 2 24 1.0e2 -t 1.0e-2 -s sp_b318_HITRAN_a16 +p -lk -o co2he_o -m co2he_m -L co2he_lbl.nc -np 16
prep_spec <<EOF
sp_b318_HITRAN_a16
a
19
y
co2he_o
-1
EOF
sed -i -e 's/\*\*\*\*\*/0/g' sp_b318_HITRAN_a16
sed -i -e 's/\*\*\*//g' sp_b318_HITRAN_a16
sed -i -e 's/            NaN/0.000000000E+00/g' sp_b318_HITRAN_a16
sed -i -e 's/          51/           0/g' sp_b318_HITRAN_a16
touch temp/done_co2he_cia
rsync -av sp_b318_HITRAN_a16 temp/sp_b318_HITRAN_a16_23
rsync -av sp_b318_HITRAN_a16_k temp/sp_b318_HITRAN_a16_23_k
echo '###############   CH4-He CIA   ###############'
Ccorr_k -F pt_grid_cia.dat -CIA /network/group/aopp/planetary/RTP015_LICHTENBERG_20GEDGE/spectral_files/dat_continua/CH4-He_2018.cia -R 1 318 -i 0.1 -ct 6 24 1.0e2 -t 1.0e-2 -s sp_b318_HITRAN_a16 +p -lk -o ch4he_o -m ch4he_m -L ch4he_lbl.nc -np 16
prep_spec <<EOF
sp_b318_HITRAN_a16
a
19
y
ch4he_o
-1
EOF
sed -i -e 's/\*\*\*\*\*/0/g' sp_b318_HITRAN_a16
sed -i -e 's/\*\*\*//g' sp_b318_HITRAN_a16
sed -i -e 's/            NaN/0.000000000E+00/g' sp_b318_HITRAN_a16
sed -i -e 's/          51/           0/g' sp_b318_HITRAN_a16
touch temp/done_ch4he_cia
rsync -av sp_b318_HITRAN_a16 temp/sp_b318_HITRAN_a16_24
rsync -av sp_b318_HITRAN_a16_k temp/sp_b318_HITRAN_a16_24_k
echo '###############   O2-CO2 CIA   ###############'
Ccorr_k -F pt_grid_cia.dat -CIA /network/group/aopp/planetary/RTP015_LICHTENBERG_20GEDGE/spectral_files/dat_continua/O2-CO2_2011.cia -R 1 318 -i 0.1 -ct 7 2 1.0e2 -t 1.0e-2 -s sp_b318_HITRAN_a16 +p -lk -o o2co2_o -m o2co2_m -L o2co2_lbl.nc -np 16
prep_spec <<EOF
sp_b318_HITRAN_a16
a
19
y
o2co2_o
-1
EOF
sed -i -e 's/\*\*\*\*\*/0/g' sp_b318_HITRAN_a16
sed -i -e 's/\*\*\*//g' sp_b318_HITRAN_a16
sed -i -e 's/            NaN/0.000000000E+00/g' sp_b318_HITRAN_a16
touch temp/done_o2co2_cia
rsync -av sp_b318_HITRAN_a16 temp/sp_b318_HITRAN_a16_25
rsync -av sp_b318_HITRAN_a16_k temp/sp_b318_HITRAN_a16_25_k
echo '###############   O2-N2 CIA   ###############'
Ccorr_k -F pt_grid_cia.dat -CIA /network/group/aopp/planetary/RTP015_LICHTENBERG_20GEDGE/spectral_files/dat_continua/O2-N2_2018_fixed.cia -R 1 318 -i 0.1 -ct 7 13 1.0e2 -t 1.0e-2 -s sp_b318_HITRAN_a16 +p -lk -o o2n2_o -m o2n2_m -L o2n2_lbl.nc -np 16
prep_spec <<EOF
sp_b318_HITRAN_a16
a
19
y
o2n2_o
-1
EOF
sed -i -e 's/\*\*\*\*\*/0/g' sp_b318_HITRAN_a16
sed -i -e 's/\*\*\*//g' sp_b318_HITRAN_a16
sed -i -e 's/            NaN/0.000000000E+00/g' sp_b318_HITRAN_a16
touch temp/done_o2n2_cia
rsync -av sp_b318_HITRAN_a16 temp/sp_b318_HITRAN_a16_26
rsync -av sp_b318_HITRAN_a16_k temp/sp_b318_HITRAN_a16_26_k
echo '###############   N2-H2O CIA   ###############'
Ccorr_k -F pt_grid_cia.dat -CIA /network/group/aopp/planetary/RTP015_LICHTENBERG_20GEDGE/spectral_files/dat_continua/N2-H2O_2018.cia -R 1 318 -i 0.1 -ct 13 1 1.0e2 -t 1.0e-2 -s sp_b318_HITRAN_a16 +p -lk -o n2h2o_o -m n2h2o_m -L n2h2o_lbl.nc -np 16
sp_b318_HITRAN_a16
a
19
y
n2h2o_o
-1
EOF
sed -i -e 's/\*\*\*\*\*/0/g' sp_b318_HITRAN_a16
sed -i -e 's/\*\*\*//g' sp_b318_HITRAN_a16
sed -i -e 's/            NaN/0.000000000E+00/g' sp_b318_HITRAN_a16
touch temp/done_n2h2o_cia
rsync -av sp_b318_HITRAN_a16 temp/sp_b318_HITRAN_a16_27
rsync -av sp_b318_HITRAN_a16_k temp/sp_b318_HITRAN_a16_27_k
echo '###############   N2-CH4 CIA   ###############'
Ccorr_k -F pt_grid_cia.dat -CIA /network/group/aopp/planetary/RTP015_LICHTENBERG_20GEDGE/spectral_files/dat_continua/N2-CH4_2011.cia -R 1 318 -i 0.1 -ct 13 6 1.0e2 -t 1.0e-2 -s sp_b318_HITRAN_a16 +p -lk -o n2ch4_o -m n2ch4_m -L n2ch4_lbl.nc -np 16
prep_spec <<EOF
sp_b318_HITRAN_a16
a
19
y
n2ch4_o
-1
EOF
sed -i -e 's/\*\*\*\*\*/0/g' sp_b318_HITRAN_a16
sed -i -e 's/\*\*\*//g' sp_b318_HITRAN_a16
sed -i -e 's/            NaN/0.000000000E+00/g' sp_b318_HITRAN_a16
touch temp/done_n2ch4_cia
rsync -av sp_b318_HITRAN_a16 temp/sp_b318_HITRAN_a16_28
rsync -av sp_b318_HITRAN_a16_k temp/sp_b318_HITRAN_a16_28_k
echo '###############   N2-H2 CIA   ###############'
Ccorr_k -F pt_grid_cia.dat -CIA /network/group/aopp/planetary/RTP015_LICHTENBERG_20GEDGE/spectral_files/dat_continua/N2-H2_2011.cia -R 1 318 -i 0.1 -ct 13 23 1.0e2 -t 1.0e-2 -s sp_b318_HITRAN_a16 +p -lk -o n2h2_o -m n2h2_m -L n2h2_lbl.nc -np 16
prep_spec <<EOF
sp_b318_HITRAN_a16
a
19
y
n2h2_o
-1
EOF
sed -i -e 's/\*\*\*\*\*/0/g' sp_b318_HITRAN_a16
sed -i -e 's/\*\*\*//g' sp_b318_HITRAN_a16
sed -i -e 's/            NaN/0.000000000E+00/g' sp_b318_HITRAN_a16
touch temp/done_n2h2_cia
rsync -av sp_b318_HITRAN_a16 temp/sp_b318_HITRAN_a16_29
rsync -av sp_b318_HITRAN_a16_k temp/sp_b318_HITRAN_a16_29_k
echo '###############   N2-He CIA   ###############'
Ccorr_k -F pt_grid_cia.dat -CIA /network/group/aopp/planetary/RTP015_LICHTENBERG_20GEDGE/spectral_files/dat_continua/N2-He_2018.cia -R 1 318 -i 0.1 -ct 13 24 1.0e2 -t 1.0e-2 -s sp_b318_HITRAN_a16 +p -lk -o n2he_o -m n2he_m -L n2he_lbl.nc -np 16
prep_spec <<EOF
sp_b318_HITRAN_a16
a
19
y
n2he_o
-1
EOF
sed -i -e 's/\*\*\*\*\*/0/g' sp_b318_HITRAN_a16
sed -i -e 's/\*\*\*//g' sp_b318_HITRAN_a16
sed -i -e 's/            NaN/0.000000000E+00/g' sp_b318_HITRAN_a16
touch temp/done_n2he_cia
rsync -av sp_b318_HITRAN_a16 temp/sp_b318_HITRAN_a16_30
rsync -av sp_b318_HITRAN_a16_k temp/sp_b318_HITRAN_a16_30_k
echo '###############   H2-CH4 CIA   ###############'
Ccorr_k -F pt_grid_cia.dat -CIA /network/group/aopp/planetary/RTP015_LICHTENBERG_20GEDGE/spectral_files/dat_continua/H2-CH4_eq_2011.cia -R 1 318 -i 0.1 -ct 23 6 1.0e2 -t 1.0e-2 -s sp_b318_HITRAN_a16 +p -lk -o h2ch4_o -m h2ch4_m -L h2ch4_lbl.nc -np 16
prep_spec <<EOF
sp_b318_HITRAN_a16
a
19
y
h2ch4_o
-1
EOF
sed -i -e 's/\*\*\*\*\*/0/g' sp_b318_HITRAN_a16
sed -i -e 's/\*\*\*//g' sp_b318_HITRAN_a16
sed -i -e 's/            NaN/0.000000000E+00/g' sp_b318_HITRAN_a16
touch temp/done_h2ch4_cia
rsync -av sp_b318_HITRAN_a16 temp/sp_b318_HITRAN_a16_31
rsync -av sp_b318_HITRAN_a16_k temp/sp_b318_HITRAN_a16_31_k
echo '###############   H2-He CIA   ###############'
Ccorr_k -F pt_grid_cia.dat -CIA /network/group/aopp/planetary/RTP015_LICHTENBERG_20GEDGE/spectral_files/dat_continua/H2-He_2011.cia -R 1 318 -i 0.1 -ct 23 24 1.0e2 -t 1.0e-2 -s sp_b318_HITRAN_a16 +p -lk -o h2he_o -m h2he_m -L h2he_lbl.nc -np 16
prep_spec <<EOF
sp_b318_HITRAN_a16
a
19
y
h2he_o
-1
EOF
sed -i -e 's/\*\*\*\*\*/0/g' sp_b318_HITRAN_a16
sed -i -e 's/\*\*\*//g' sp_b318_HITRAN_a16
sed -i -e 's/            NaN/0.000000000E+00/g' sp_b318_HITRAN_a16
touch temp/done_h2he_cia
rsync -av sp_b318_HITRAN_a16 temp/sp_b318_HITRAN_a16_32
rsync -av sp_b318_HITRAN_a16_k temp/sp_b318_HITRAN_a16_32_k
rsync -av sp_b318_HITRAN_a16 temp/sp_b318_HITRAN_a16_no_spectrum
rsync -av sp_b318_HITRAN_a16_k temp/sp_b318_HITRAN_a16_no_spectrum_k
prep_spec <<EOF
sp_b318_HITRAN_a16
a
6
n
T
100 4000
100
2
n
/home/lichtenberg/codes/socrates/socrates_2002/data/solar/kurucz_95
y
-1
EOF
