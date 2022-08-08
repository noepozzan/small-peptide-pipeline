#!/bin/bash

# for now, instead use this direct python call by running it in the ${projectDir} directory
# ATTENTION: needs >=python3.7 to run
# ALSO: might have to adapt the config files:
# I had to add test profiles bc gitlab runner doesn't have slurm installed

mkdir -p ~/.singularity/cache/library

python ${PWD}/data/scripts/pull_containers.py \
	--slurm_in ${PWD}/conf/envs/slurm.config \
	--slurm_out ${PWD}/conf/envs/slurm_offline.config \
    --singularity_in ${PWD}/conf/envs/singularity.config \
    --singularity_out ${PWD}/conf/envs/singularity_offline.config \
	--dest ${HOME}/.singularity/cache/library/

