set -e
# 1.) source this file instead of executing it
#     source install_developer_conda.sh
#
# This script is just setting up a development environment
# (in the sense of setup.py develop)
# It is NOT trying to replace a proper conda package.
# Since we frequently will have to work on bgc_md2 and
# CompartmentalSystems, LAPM, testinfrastructure simultaneuously we also
# install them in development mode
# To make that easy we assume that there repositories have been checked
# out under src/ComartmentalSystems src/LAPM and src/testinfrastructure 

envname="bgc_md2"
conda remove --name $envname --all
conda env create  -f environment.yml
conda activate $envname

for dir in CompartmentalSystems LAPM testinfrastructure
do 
  cd src/${dir}
  source install_developer_conda.sh
  cd -
done
python setup.py develop
