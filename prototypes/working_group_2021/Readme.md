# Unified data assimilation resources for the working group
## How to get the code and contribute:
This folder is part of a bigger repository.
To use the code You have two options. (Please read both paragraph before you decide ;-) )
1. Clone the repository 
   * Pro:
      * I will grant you the priveliges to push which allows you to make changes available quicker to everybody
        without any interaction with me.
      * Since you will most likely (99% of the time) work in your model's folder the chances of collisions are 
        small.
      * If you make a muscle memory habit of pulling the changes before you start your own work conflicts are
        even less likely
      * If you commit your changes to a temporary branch like `test` or one of your own the master or main 
        branch of the repository (which is used by all the other users) is still save.
      


    * Con:
      If you don't know what a merge conflict is or how to resolve it without messing up somebody else's work or tend to push to the `master` branch without thinking (and testing ;-) before or if you are really really scared then you are better off creating your own fork...
      with you requesting us to pull your changes...
    * How to clone: 
      1. Go to the main site of the repo: https://github.com/MPIBGC-TEE/bgc_md2
      2. Press the green Code button and copy the url
      3. type 
         ```bash
         git clone https://github.com/MPIBGC-TEE/bgc_md2.git --recurse submodules
         ``` 
         in your terminal 
         (in the directory where you want to checkout the code)
      
1. Fork the repository:
    * Pro:
      Your repo is more isolated:
      Forking (creating your own branch of the repo) has the advantage that it is harder to mess up anybody else's work,
      even if you are not too familiar with git, since you just cannot directly push to the original repo.
    * Con:
      Your repo is more isolated ;-):
      * To make your changes available to the group you have to make a pull request and I have to pull your code.
      * Updating your fork with changes from the main repository is slightly more tedious and requires 
        more steps than updating some branch of the main repo. (Although the process is well documented on github) 
      * If you have trouble managing your fork, you will have to invite me as a colllaborator to be able to help	you (since I have no priveliges on your fork)
      
      
    * How to fork
      1. read the general github instructions here: https://docs.github.com/en/get-started/quickstart/fork-a-repo 
      1. press the Fork button at the upper right cornder of your github window.
      1. Clone you own fork



## How to use the code
### Structure of the directory
The files are organized in two catagories:
 * model specific
 * general
For every model there is a specific sub folder, containing the model specific files, while the general functions reside one level above.

We have created one model that examplifies the workflow you will follow to add new models. 
The folder name is  `bgc_md2/prototypes/working_group_2021/kv_visit2`
At the beginning you will only need one file in this folder.
`createModel.py`

### Test the example

Before you start working on your own model you should make sure that you can run the example.
1. Install the `bgc_md2` package as described in the instructions in the gitHub https://github.com/MPIBGC-TEE/bgc_md2#installation
   Remark: If you cloned the repo change the url to your fork, otherwise you will not see the effect of your latest codechanges immediately.

1. `cd kv_visit2` 
   ```bash
   jupyter notebook
   ```
   This will fire up a `jupyter` server and most likely open a browser window.
   If it does not you can point your browser to `http://localhost:8888`
1. In the browser window open the file `createModel.py`
   If it opens as a jupyter notebook then you already have `jupytext` installed.
   It it doesn't you have to install jupytext 
   https://jupytext.readthedocs.io/en/latest/install.html

   (Ohterwise we will have to do with notebook (`*.ipynb`) files which are difficult to handle by 
   version control systems)


### Adapt the code to your own needs
1. Create a folder for your model. We call it `{your_model}` from now on.
1. Copy the `createModel.py` file into your new folder and open it with `jupyter`
1. Start changing the code as described in the notebook.
1. Later we will outsource the code into different files and streamline the notebooks a bit and add more functionality.

### General Remarks and Additional information


