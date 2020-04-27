# Where are the image folders?
IMAGE_DIR_BASE=/Users/tim/bitbucket/pcd_couple-interior-atmosphere/atm_rad_conv/output/mov_figs

# Where to save the movie files
MOVIE_DIR=/Users/tim/bitbucket/pcd_couple-interior-atmosphere/atm_rad_conv/output

# Frames per second
FPS=24

# Movie length [s]
ML=60

# Loop through simulation branches
# BATCH=test20_moonmars_pebble2 # test15_mars025 test16_diapirs test19_moonmars_pebble1 test20_moonmars_pebble2
# test21_moonmars_imp1 test22_moonmars_imp2 test23_moonmars_pebble3 test24_moonmars_imp3
for BATCH in trpp_H2O; do

    # Go the base directory of the folders which contain the images
    cd ${IMAGE_DIR_BASE}/${BATCH}
    pwd

    #for DIR in $(ls -d r*/); do
    for FIELD in $BATCH; do

        # Get movie name from folder name of the images
        # MOVIE_NAME=${DIR%%/}
        MOVIE_NAME=${PWD##*/}

        # Absolute directory of where images are stored
        # IMAGE_DIR=${IMAGE_DIR_BASE}/${BATCH}/${MOVIE_NAME}
        IMAGE_DIR=${IMAGE_DIR_BASE}/${BATCH}
        # MOVIE_NAME=${PWD##*/}

        # Calculate shown images per second, such that movies are always $ML seconds long
        IMAGE_NO=$(ls $IMAGE_DIR | wc -l)
        IPS=$(echo "$IMAGE_NO/$ML" | bc)

        # Inform user about the input files
        echo "Use the following files as input:"
        ls $IMAGE_DIR/*${FIELD}*.png

        # Process the files and create movie
        ffmpeg -framerate $IPS -pattern_type glob -i "$IMAGE_DIR/${FIELD}*.png" -vf scale=1504:1232 -y -c:v libx264 -r $FPS -pix_fmt yuv420p -profile:v baseline -level 3.0 -movflags +faststart -preset veryslow -crf 18 ${MOVIE_DIR}/${MOVIE_NAME}_${FIELD}.mp4

    # Lossless video: https://trac.ffmpeg.org/wiki/Encode/H.264
    done
done
