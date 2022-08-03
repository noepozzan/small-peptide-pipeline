#!/bin/bash

# first of all, set some env variables
# the environment variables will be needed to access the private docker registry
cat ./data/bash_scripts/echo_env.sh
source ./data/bash_scripts/echo_env.sh

# convert proteomics raw to mzML files
nextflow run \
	main.nf \
	-profile slurm \
	-entry RAW_TO_MZML

# run this to pull images from image repos to local
# this ensures error-free execution of the pipeline (atleast at this step, hopefully)
bash ./data/scripts/pull_containers.sh

# we might also want to pull the newest fasta genome and annotation files:
# This subworkflow takes a special profile, since it has to run on the login node
# which is the only one connected to the internet
# With this little workflow, it should probably be possible to download the files that you would like.
# But honestly implementing this nicely would probably take some time, so I'll think about it a later point.
nextflow run \
	main.nf \
	-profile full,pull

nextflow run main.nf \
	-profile full,cluster

# generate test files by running
: '
for riboseq fastq files:
zcat ${riboseq_file} \                                                      
        | head -n ${lines} \                                                    
        | gzip -c > TEST_${riboseq_file}

# after lots of trying and tempering with many tools to get a small mzML file that can pass the integration test
# (a file that is not too small), I have found MzMLSplitter, a tool available through msopen v2
# singularity shell docker://blcdsdockerregistry/openms-thirdparty:2.6.0
# cd /path/to/mzML files
# mkdir parts
# MzMLSplitter -in 010_TSP7-21_D19.mzML -parts 10 -out parts/010

the rest of the files (fasta files etc..) can just be cut with "head"
'

# if your proteomics filenames do not clearly indicate the experimental run
# and are not in the right format, this process will fix it.
# In order for this to work, the file experimental_conditions.csv in the data folder
# has to be quickly filled out.
nextflow run \
	main.nf \
	-entry MAP_NAMES \
	-profile slurm_offline

# integration test

nextflow run \
    main.nf \
    -profile test,slurm

# finally call sbatch to make the full main workflow run
sbatch slurm_script
# or just: nextflow run main.nf -profile full,slurm_offline
