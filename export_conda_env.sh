# this file wants to be sourced! (otherwise the conda activate does not work)
conda config --add channels anaconda
conda config --add channels conda-forge
conda update -y -n base -c defaults conda
conda activate base
conda env remove -n binder
conda env create -f environment_free.yml && conda activate binder && ./postBuild && conda env export > environment.yml
