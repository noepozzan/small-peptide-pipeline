#!/bin/bash

# get dir to correctly change dir later
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

# make sure git lfs works
cd ~
git lfs install --skip-repo

# change into the project's main dir
cd $SCRIPT_DIR/../..

# get the test data and clean after it
git clone https://github.com/noepozzan/small_peptide_pipeline_test_data.git
mkdir -p data/tests/
mv small_peptide_pipeline_test_data/* data/tests/
rm -rf small_peptide_pipeline_test_data

