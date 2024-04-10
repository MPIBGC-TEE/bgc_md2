# automatically  procuded by scripts/make_conda_clone_installers.py , do not change manually
micromamba install -y -c conda-forge --file requirements.test --file requirements.doc --file requirements.non_src pip
pip install -r requirements.src
