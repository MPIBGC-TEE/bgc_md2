
## Installation
We assume here that you have a recent installation of python3 and conda. 
For developers who work with CompartmentalSystems LAPM and testinfrastructure simultaneously: 
   * Clone the repository and its [submodules](https://git-scm.com/book/en/v2/Git-Tools-Submodules):
   ```
   git clone --recurse-submodules https://github.com/MPIBGC-TEE/bgc_md2.git 
   ```
   * To pull the changes in bgc_md2 and the submodules simultaneuously use:
   ```
   git pull --recurse-submodules
   ```
   (Or configure an alias ```spull``` by)
   ```
   git config alias.spull 'pull --recurse-submodules'
   ```
   * Make sure that the submodule folders in `src` are not empty. In case they are you need to run
   ```
   git submodule init
   git submodule update
   ```
   * create a conda environment and run the install script    
   ```bash 
   conda create -y --name bgc_md2
   conda activate bgc_md2
   ./install_developer_conda.sh 
   ```
   This will install the dependencies and run ```python setup.py develop``` for every subpackage so that your code changes  in one of these packages take immediate effect.
   
## Structure
This prototype represents models as subfolders of folder 
```bash
models
```
## Documentation
* The latest build of the package documentation can be found [here:](https://mpibgc-tee.github.io/bgc_md2/).


## Objectives
The package is supposed to assist in the creation of 'reports' in the form of jupyter notebooks.
The notebooks will be of the following types.
1. Investigations of a single model (or modelrun).
1. Comparisons between models (modelruns).

In the first case the role of the `bgc_md` package is to guide the user (=author of a particular notebook concerned with a particular model, and simultaniously author of the `source.py` of that model) by using the computability graph (as represented by `bgc_md/resolve/MvarsAndComputers.py`) to either
* show which addidional results can be computed, given the information already present in the models `source.py' or
* show which additional information has to be provided in the models `source.py` to be able to obtain a desired result.

In the second case the same assistance is required for queries, which are best described by examples. 
* Create a table including all the models for which we can compute the compartmental matrix (from whatever Mvars are provided, in the different model files)
* Compute the maximum set of `Mvars` we can compute for a given set of models
* ...



## various notes on implementation

* The 'Computers' and 'MVars' represent a set of types and strictly typed
  functions (including the return values).
  This has been implemented with the new python type annotations.  
  An advantage is that we can express our
  idea in a well defined and well documented way and avoid extra effort for the
  user..  
* The computibility graph is expensive to create and only changes if new
  `Computers` and `MVars` are created.  It should be cached, which encurages
  the use of immutable data structures. (since we can use functools )

   


