# automatically  procuded by scripts/make_conda_clone_installers.py , do not change manually
#
# We install as much as we can via mamba and leave only the src packages for pip.
mamba install -y -c conda-forge --file requirements.test --file requirements.doc --file requirements.non_src pip
pip install -r requirements.only_pip
pip install -r requirements.src
