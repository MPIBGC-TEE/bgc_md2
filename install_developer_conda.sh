#!/bin/bash
set -e
# The script is supposed to run in an activated conda environment.
# e.g. after 
# conda create -y --name bgc_md2 
# conda activate bgc_md2
# and in this directory
#
# It is NOT trying to replace a proper conda package.
# Since we frequently will have to work on bgc_md2 and
# CompartmentalSystems, LAPM, testinfrastructure simultaneuously we also
# install them in development mode
# To make that easy we assume that there repositories have been checked
# out under src/ComartmentalSystems src/LAPM and src/testinfrastructure 

conda install -y -c conda-forge python=3.10 #should install the tested version of python3
for dir in testinfrastructure ComputabilityGraphs CompartmentalSystems LAPM 
do 
  echo '#################'
  echo $dir
  cd src/${dir}
  ./install_developer_conda.sh
  cd -
done
conda install -y -c conda-forge --file requirements.conda
python setup.py develop
