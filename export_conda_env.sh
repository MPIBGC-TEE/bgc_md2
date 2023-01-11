conda env remove -n binder
conda env create -f environment.yml && conda activate binder && ./postBuild && conda env export > environment.yml
