# refineSetFields.sh
#
# setFields on a refined grid (which will be what the interface "sees" when
# running, which gives a more accurate initial condition
#
# Requires system/setFieldsDictInit which sets the initial area to refine over
#
# Author:   Yatin Darbar (CDT in Fluid Dynamics, University of Leeds, UK)
# Made for: OpenFOAM Foundation 9
# Updated:  10 August 2022
# Version:  2.0

echo -e "refineSetFields (version 2.0)\nMade for OpenFOAM 9\n"

# Check the setFieldsDictInit dictionary is accessible
if [ ! -f system/setFieldsDictInit ]; then
    echo "The initial set fields domain should be defined in" \
         "system/setFieldsDictInit. This file is currently not" \
         "accessible - check it!"
    exit 1
fi

# Create log file to redirect terminal output
logName=rsfLog$(date -d " today" +"_%Y%m%d-%H.%M")
exec &> $logName

# Set up folder/file structure to create refined mesh
folderName=RSF$(date -d " today" +"_%Y%m%d-%H.%M")
mkdir $folderName
cp -r 0/ constant/ system/ $folderName/
cd $folderName

# Replace the usual setFields dictionary with the initial version
mv system/setFieldsDictInit system/setFieldsDict

# Get number of refinements (each direction refined by 2^maxRefine)
maxRefine=$(grep -E 'maxRefinement.*;' constant/dynamicMeshDict | \
            sed 's/;.*//' | \
            sed 's/^.*maxRefinement//' | \
            sed 's/^.* //' | \
	    sort -rn | \
	    head -1 \
           )
# Increment maxRefine by 1 to be sure (maxRefine = maxRefine + 1)
let maxRefine++    

# Set up controlDict for generating refined mesh
# Protect "stopAt endTime;" against the sed for endTime
sed -i '/stopAt/{s/.*/stopAt            placeholder;/}' system/controlDict
# Define endTime such that enough iterations are carried out to fully refine
sed -i "/endTime/{s/.*/endTime           ${maxRefine}e-7;/}" system/controlDict
# Replace/confirm endTime
sed -i '/stopAt/{s/.*/stopAt            endTime;/}' system/controlDict
# Control the number of iterations, with respect to the endTime
sed -i '/deltaT/{s/.*/deltaT            1e-7;/}' system/controlDict
# Write at (known) time steps so only at the end of the simulation
sed -i '/writeControl/{s/.*/writeControl      timeStep;/}' system/controlDict
# Write only at the end when all refining is complete (saves time)
sed -i '/writeInterval/{s/.*/writeInterval     1;/}' system/controlDict
# Don't compress so polyMesh is readily accessible 
sed -i '/writeCompression/{s/.*/writeCompression  uncompressed;/}' system/controlDict
# Ensure time step is not changed so number of iterations is known
sed -i '/adjustTimeStep/{s/.*/adjustTimeStep    no;/}' system/controlDict

# Minimises the number of iterations required to fully refine the mesh
sed -i '/refineInterval/{s/.*/refineInterval  1;/}' constant/dynamicMeshDict

# Run commands for to generate the refined initial mesh
# blockMesh is the same for the original case, so save the output separately
blockMesh > bLog
setFields
diffusiveInterFoam -noFunctionObjects

# Find (last) output time (by last file/directory written)
time=$(ls -t | head -n 1)

# Get mesh into original case 
cp -r $time/polyMesh ../constant

# Clean up
mv bLog ..
cd ..
rm -r $folderName

# Set Fields on the refined mesh
setFields > sLog
