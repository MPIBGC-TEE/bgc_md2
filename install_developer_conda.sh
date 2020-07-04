# source this file instead of executing it
# source install_developer_conda.sh
conda remove --name bgc_md2_env --all
conda env create -f environment.yml
conda activate bgc_md2_env
pip install -e .
