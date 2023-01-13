mamba update
mamba env remove -n binder
mamba env create -f environment_free.yml && conda activate binder && ./postBuild && mamba env export > environment.yml
