#!/bin/bash

# for now, instead use this direct python call by running it in the ${projectDir} directory
# ATTENTION: needs >=python3.7 to run
# ALSO: might have to adapt the config files:
# I had to add test profiles bc gitlab runner doesn't have slurm installed

python ${PWD}/data/scripts/pull_containers.py \
	--config ${PWD}/conf/slurm.config \
	--out ${PWD}/conf/slurm_offline.config \
	--dest ${HOME}/.singularity/cache/library/

