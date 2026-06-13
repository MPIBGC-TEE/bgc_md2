# This script builds a list of packages for a conda install command
# This is useful if you want to use bgc_md2 alongside other packages in a conda environment.
# It uses the full environment.yml file created for binder
# which contains two sources of package name
# 1.) the list of packages installed by conda
# 2.) the list of packages handled by pip (which evaluates the "setup.py" of our source packages to find their dependencies)

# We want a file that can be used 
import yaml 
from pathlib import Path
with Path("environment.yml").open("r") as stream:
    yd = yaml.safe_load(stream)
# find the section in the environment that specifies dependencies handled by pip
# Assuming that we are in a conda environment we don't want those
# dependencies actually handled by pip but by conda
deplist = yd["dependencies"]
conda_strings = deplist[:-1] #list element is a dict
pip_strings = deplist[-1]["pip"] #list element is a dict

names = [s.split("=")[0] for s in conda_strings] + [s.split("==")[0] for s in pip_strings]   

with Path("requirements.mm").open("w") as f: 
    f.write("\n".join(names))

