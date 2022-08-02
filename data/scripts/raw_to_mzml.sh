#!/bin/bash

# the idea for this procedure comes from:
# https://sourceforge.net/p/proteowizard/mailman/proteowizard-support/?limit=250&viewmonth=201907
# then search for "singularity"

# set some variables
files="$@"
#delete_mode="$2"
workd=$(pwd)
# create a temp dir and move into it
rm -rf tmp_dir
mkdir tmp_dir
for file in $files
do
	cp $file tmp_dir
done
cd tmp_dir
# run singularity (must be available), maybe check beforehand 
singularity build --sandbox pwiz_box docker://chambm/pwiz-skyline-i-agree-to-the-vendor-licenses
# wineprefix64 has to be moved to user's home due to ownership issues..
mv pwiz_box/wineprefix64/ ~
# regenerate an empty wineprefix64 inside the singularity file
mkdir pwiz_box/wineprefix64
# build the singularity image again
singularity build singularity_pwiz.sif pwiz_box
# run the singularity container
for file in $files
do
	singularity exec -B ~/wineprefix64:/wineprefix64 ./singularity_pwiz.sif wine msconvert $file
done
# move mzML file out of temp dir
mv *.mzML $workd
# function to remove temp dir
cleanup () {
    rc=$?
    cd $workd
	rm -rf tmp_dir
    echo "Exit status: $rc"
}
if [[ "$delete_mode" == "true" ]]
then
	cleanup
fi

