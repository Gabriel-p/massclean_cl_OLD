#!/bin/bash


##  MASSCLEAN V2.0123 ##

# This script generates the synthetic clusters for the specified initial mass distance,
# age and metallicity ranges and stores the files in the corresponding folders.

# The first run with the smallest initial mass should be done manually to also generate
# the field stars field.plot files for each distance and A_V set of parameters.



# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# 1- The 'massclean2.013' folder should be a sub-folder of the one where this
#    script is located.
# 2- The convention for the names of the created folders is xx_yyyy where 'xx' is the initial mass
#    of the clusters (2 chars) and 'yyyy' is the distance in parsecs (3 or 4 chars)



# Current directory.
CDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
# Location of the cluster.ini file.
CLUSTER_INI=$CDIR"/massclean2.013/ini.files/cluster.ini"

# Change to code directory.
cd $CDIR"/massclean2.013"

# Clean left over files from previous run.
./clean.all

# To generate field.plot files. Only run these commands once (for the smallest
# initial mass so it will be faster) since the rest of the clusters will use
# these same fields.
#./clean.all
#./all.run
#./goimage2
#./gofield2
#chmod u+x add_field
#./add_field


# Declare arrays of metallicities, ages, initial masses and distance.
METAL=('002' '008' '019' '030')
AGES=('0600' '0750' '0900' '1000')
INIT_MASS=('02' '04' '06' '08')
INIT_MASS2=('200' '400' '600' '800')
DIST=('500' '1000' '1500' '2000')
AV=('0.5' '1.0' '1.5' '2.0')
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
            sed -i "45s/.*/${INIT_MASS2[${j}]}    (2)/" $CLUSTER_INI
            sed -i "46s/.*/${DIST[${k}]}    (3)/" $CLUSTER_INI
            sed -i "56s/.*/${AV[${k}]}    (13)/" $CLUSTER_INI

            #read -p "Press [Enter] key to continue..."
            
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
            CODE_DIR="./clusters/${INIT_MASS[${j}]}_${DIST[${k}]}"
            # Create sub-folder if it does not exist.
            mkdir -p $CODE_DIR

            #read -p "1 Press [Enter] key to continue..."

            # Copy synth clusters files to output folder.
            cp ./data/iso.trek/*.trek $CODE_DIR
            # Copy field stars file to output folder.
            cp ./data/field.plot $CODE_DIR
            # Copy used cluster.ini file to folder.
            cp $CLUSTER_INI $CODE_DIR

        done
    done
done

echo 'Done.'
