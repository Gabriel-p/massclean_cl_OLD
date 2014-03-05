#!/bin/bash



# This script generates the synthetic clusters for the specified initial mass and distance ranges
# and stores the files in the corresponding folders. The first run with the smallest initial
# mass should be done manually to also generate the field stars field.plot files for each
# distance and A_V set of parameters.



# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# 1- This script should be located inside the 'massclean2.012' folder.
# 2- The sub folders inside the CODE_DIR directory should already be created since
# the script will not work otherwise. The convention is xx_yyyy where 'xx' is the initial mass
# of the clusters (2 chars) and 'yyyy' is the distance in parsecs (3 or 4 chars)




# To generate field.plot files. Only run these commands once (for the smallest
# initial mass so it will be faster) since the rest of the clusters will use
# these same fields.
#./clean.all
#./all.run
#./goimage2
#./gofield2
#chmod u+x add_field
#./add_field


# Location of the cluster.ini file.
CLUSTER_INI=~/Descargas/massclean2.012/ini.files/cluster.ini

#Declare arrays of initial mass and distance.
INIT_MASS=('01' '02' '03' '04' '05' '06' '07' '08')
INIT_MASS2=('100' '200' '300' '400' '500' '600' '700' '800')
DIST=('500' '1000' '1500' '2000' '2500' '3000')
AV=('0.5' '1.0' '1.5' '2.0' '2.5' '3.0')
# get number of elements in the initial mass array
MASSES=${#INIT_MASS[@]}
# Same for distances.
DISTANCES=${#DIST[@]}

# Iterate through all initial masses.
for (( i=0;i<$MASSES;i++)); do
    
    # Iterate through all distances.
    for (( j=0;j<$DISTANCES;j++)); do
        #echo ${INIT_MASS[${i}]}, ${DIST[${j}]}
        
        # Modify 'cluster.ini' file.
        sed -i "40s/.*/${INIT_MASS2[${i}]}    (2)/" $CLUSTER_INI
        sed -i "41s/.*/${DIST[${j}]}    (3)/" $CLUSTER_INI
        sed -i "51s/.*/${AV[${j}]}    (13)/" $CLUSTER_INI
        
        # Generate clusters with different ages according to the parameters set above.
        ./clean.all
        ./all.run
        
        # Set dir to move the files to.
        CODE_DIR="/media/rest/Dropbox/GABRIEL/CARRERA/3-POS-DOC/trabajo/codigo/massclean/${INIT_MASS[${i}]}_${DIST[${j}]}"
        echo "${INIT_MASS[${i}]}_${DIST[${j}]}"
        # Define cluster's ages to copy.
        file1=~/Descargas/massclean2.012/data/iso.trek/is1_p019_0650.trek
        file2=~/Descargas/massclean2.012/data/iso.trek/is1_p019_0700.trek
        file3=~/Descargas/massclean2.012/data/iso.trek/is1_p019_0750.trek
        file4=~/Descargas/massclean2.012/data/iso.trek/is1_p019_0800.trek
        file5=~/Descargas/massclean2.012/data/iso.trek/is1_p019_0850.trek
        file6=~/Descargas/massclean2.012/data/iso.trek/is1_p019_0900.trek
        file7=~/Descargas/massclean2.012/data/iso.trek/is1_p019_0950.trek
        # Move cluster files to code folder.
        cp $file1 $file2 $file3 $file4 $file5 $file6 $file7 $CODE_DIR
        # Copy used cluster.ini file to folder.
        cp $CLUSTER_INI $CODE_DIR
        
    done 
    
done 

echo 'Done'
