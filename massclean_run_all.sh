#!/bin/bash


##  MASSCLEAN V2.0123 ##

# This script generates the synthetic clusters for the specified initial mass,
# distance, age and metallicity ranges and stores the files in the
# corresponding folders.


# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# 1- The 'massclean2.01X' folder should be a sub-folder of the one where this
#    script is located.
# 2- The convention for the names of the created folders is xx_yyyy where 'xx'
#    is the initial mass of the clusters (2 chars) and 'yyyy' is the distance in
#    parsecs (3 or 4 chars)
# 3- A field plot is generated for each distance although this could probably
#    be improved by moving the ./gofield2 command outside of the for loops.
# 4- The file 'field.ini' should use the J band as reference in (1) since
#    apparently the V band causes troubles. The max_J mag can be set as:
#    max_J = max_V - 3 in (5). 



# Current directory.
CDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
# MASSCLEAN version
ver="2.014"
# Location of the cluster.ini file.
CLUSTER_INI=$CDIR"/massclean${ver}/ini.files/cluster.ini"
KING_INI=$CDIR"/massclean${ver}/ini.files/king.ini"
FIELD_INI=$CDIR"/massclean${ver}/ini.files/field.ini"

# Change to code directory.
cd $CDIR"/massclean${ver}"

# Remove *ALL* previous synthetic clusters created.
rm -rfv clusters/*

# Clean left over files from previous run.
./clean.all

# Set parameters for field stars.
# Use J band as per Popescus' recommendation.
sed -i "34s/.*/J         (1)/" $FIELD_INI
# Number of field stars.
sed -i "36s/.*/12000      (3)/" $FIELD_INI
# Minimum magnitude.
sed -i "37s/.*/22        (4)/" $FIELD_INI
# Maximum magnitude.
sed -i "38s/.*/8         (5)/" $FIELD_INI

# Set field size.
sed -i "47s/.*/4.8    (4)/" $CLUSTER_INI
sed -i "50s/.*/2048   (7)/" $CLUSTER_INI
sed -i "51s/.*/2048   (8)/" $CLUSTER_INI

# Set cluster r_t to 250px: 250/(2048/2) = 0.244
sed -i "36s/.*/0.244    (3)/" $KING_INI

# Declare arrays of metallicities, ages, initial masses and distance.
# METAL=('002' '008' '019' '030')
# AGES=('0700' '0800' '0900')
# INIT_MASS=('50' '100' '200' '300' '400' '500' '600' '800' '1000')
# DIST=('500' '1000' '3000' '5000')
# AV=('0.1' '0.5' '1.0' '3.0')

METAL=('015')
AGES=('0800')
INIT_MASS=('1000')
DIST=('5000')
AV=('0.5')

# Get number of elements in the initial metallicity array.
METAL_n=${#METAL[@]}
# Same for ages.
AGES_n=${#AGES[@]}
# Same for masses.
MASSES_n=${#INIT_MASS[@]}
# Same for distances.
DISTANCES_n=${#DIST[@]}

# Iterate through all initial masses.
for (( i=0;i<$METAL_n;i++)); do

    # Use set.ini file for standard UBVI photometry.
    cp ./ini.files/default/sets.ini1  ./ini.files/sets.ini
    # Remove previously used isochrones.
    rm ./isochrones/iso.to.use/*.*
    # Select which isochrones to use (metallicity and age)
    for (( h=0;h<$AGES_n;h++)); do
        # Define isoch to use.
        file_a="./isochrones/padova/${METAL[${i}]}/is1_p${METAL[${i}]}_${AGES[${h}]}.PADOVA"
        # Copy it to /iso.to.use folder.
        cp $file_a ./isochrones/iso.to.use
    done

    # Iterate through all initial masses.
    for (( j=0;j<$MASSES_n;j++)); do
        
        # Iterate through all distances.
        for (( k=0;k<$DISTANCES_n;k++)); do
            #echo ${METAL[${i}]}, ${INIT_MASS[${j}]}, ${DIST[${k}]}
            
            # Modify 'cluster.ini' file.
            sed -i "45s/.*/${INIT_MASS[${j}]}    (2)/" $CLUSTER_INI
            sed -i "46s/.*/${DIST[${k}]}    (3)/" $CLUSTER_INI
            sed -i "56s/.*/${AV[${k}]}    (13)/" $CLUSTER_INI
            
            # Generate clusters with different ages and metallicities
            # according to the parameters set above.

            # Clean left over files from previous run.
            ./clean.all

            # Sequence to replace ./all.run.
            ./goimf2
            ./goking2
            ./goccm2

            # Run final script to generate the synth clusters.
            # Uses all the isochrones stored in folder /iso.to.use
            ./gotrek2

            # Generate field of stars.
            ./gofield2
            
            # Set dir to move the synthetic clusters to.
            CODE_DIR="./clusters/${INIT_MASS[${j}]}_${DIST[${k}]}_${AV[${k}]}"
            # Create sub-folder if it does not exist.
            mkdir -p $CODE_DIR

            # echo $CODE_DIR
            # read -p "1 Press [Enter] key to continue..."

            # Copy synth clusters files to output folder.
            cp ./data/iso.trek/*.trek $CODE_DIR
            # Copy field stars file to output folder.
            cp ./data/field.plot $CODE_DIR
            # Copy used cluster.ini and field.ini files to folder.
            cp $CLUSTER_INI $CODE_DIR
            cp $FIELD_INI $CODE_DIR

        done
    done
done

# Change back to main directory.
cd $CDIR
# Run python script to merge the clusters and the fields.
python massclean_cluster_field_merge.py

echo 'Done.'
