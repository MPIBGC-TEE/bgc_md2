conda update -y -n base -c defaults conda
mamba env remove -n binder
mamba env create -f environment_free.yml && conda activate binder && ./postBuild && mamba env export > environment.yml
