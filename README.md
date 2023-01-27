![test_conda_install](https://github.com/MPIBGC-TEE/bgc_md2/workflows/test_conda_developer_installation.yml/badge.svg)
![test_conda_install](https://github.com/MPIBGC-TEE/bgc_md2/workflows/test_conda_binder_pinned.yml/badge.svg)
![test_conda_install](https://github.com/MPIBGC-TEE/bgc_md2/workflows/test_conda_binder_free.yml/badge.svg)

![test_debian_pip_install](https://github.com/MPIBGC-TEE/bgc_md2/workflows/test_debian_pip_install.yml/badge.svg)

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/MPIBGC-TEE/bgc_md2/binder)
## Installation

Please read carefully before you  "Copy and paste", since some instructions only make sense for certain platforms.

* For developers who work with CompartmentalSystems LAPM and testinfrastructure simultaneously: 
   * Clone the repository and its [submodules](https://git-scm.com/book/en/v2/Git-Tools-Submodules):
     Note: 
     If you have forked the repository you probably want to put your forks ulr here.
     You can of course install the original but especially if you do so using setuptools `develop` mode
     to see the effects of your changes on the installed package you probably want to install from your fork.

     * If you do *not* have a bgc_md2 repo yet:
       ```bash
       git clone --recurse-submodules https://github.com/MPIBGC-TEE/bgc_md2.git
       ```

       
     * If you *already* have a bgc_md2 repo (and want to keep it):
        * Pull the changes in bgc_md2 and the submodules simultaneuously:
          ```
          git pull --recurse-submodules
          ```
        * Make sure that the submodule folders in `src` are not empty. 
          ```
          git submodule init
          git submodule update
          ```
   * Update conda
     ```
     conda update --all
     ```
   * Create a conda environment and run the install script:
     ```bash 
     conda create -y --name bgc_md2 python=3
     conda activate bgc_md2
     cd bgc_md2
     ./install_developer_conda.sh 
     ```
     (on MS-Windows replace the last line with)
     ```bash
     install_developer_conda.bat 
     ```
     This will install the dependencies and run ```python setup.py develop``` for every subpackage so that your code changes 
     in one of these packages take immediate effect.
     
   * Run the tests.
      ```
      cd tests
      ./run_tests.py
      ```
     (on MS-Windows replace the last line with)
     ```
     python run_tests.py
     ```
      If you can run this script successfully, you have a working installation of bgc_md and can run all functions. 
  
   * Troubleshooting:
      * We noticed that in MacOS, it is necessary to update packages in the conda environment before running the tests successfully.
        Try to update conda ( ```conda update --all)``` and run the tests again.
        
   * Working with the installation:
      * pulling:
        Since you will nearly always pull with the ```--recurse-submodules``` flag   
        consider creating an alias
        ```
        git config alias.spull 'pull --recurse-submodules'
        ```
        which enables you to say  ```git spull``` to achieve the same effect
        
      * Tips to work with [git submodules:](https://git-scm.com/book/en/v2/Git-Tools-Submodules)
   

## Documentation
* The latest build of the package documentation can be found [here:](https://mpibgc-tee.github.io/bgc_md2/).


## Objectives
The package is supposed to assist in the creation of 'reports' in the form of jupyter notebooks.
The notebooks will be of the following types.
1. Investigations of a single model (or modelrun).
1. Comparisons between models/modelruns.

In the first case the role of the `bgc_md` package is to guide the user (=author of a particular notebook concerned with a particular model, and simultaniously author of the `source.py` of that model) by using the computability graph (as represented by `bgc_md/resolve/MvarsAndComputers.py`) to either
* show which addidional results can be computed, given the information already present in the models `source.py' or
* show which additional information has to be provided in the models `source.py` to be able to obtain a desired result.

In the second case the same assistance is required for queries, which are best described by examples. 
* Create a table including all the models for which we can compute the compartmental matrix (from whatever Mvars are provided, in the different model files)
* Compute the maximum set of `Mvars` we can compute for a given set of models
* ...

# Contribution
We try to keep the master green and develop new features or bug fixes in short lived branches that are then 
merged back into the master https://nvie.com/posts/a-successful-git-branching-model/
See also https://git-scm.com/book/en/v2/Git-Branching-Branches-in-a-Nutshell
https://git-scm.com/book/en/v2/Git-Tools-Submodules for the work on the dependencies.

* Example workflow to work on a feature branch `iss26-non-importable-models` you are asked to review 
  * `git branch -a` (shows all branches including remotes)
  * `git checkout --track origin/iss26-non-importable-models` (creates a local copy that you can test)
* Example to create your own feature branch (here to fix an issue )
  * `git checkout -b iss53`

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
