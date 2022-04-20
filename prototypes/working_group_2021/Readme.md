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
   
1. Update your repository everytime before you make any changes (see below)

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
1. Also open the new file `inspectModel.py` this is how the notebook looks after the refactoring.


### Adapt the code to your own needs
1. Create a folder for your model. We call it `{your_model}` from now on.
1. Copy the `createModel.py` file into your new folder and open it with `jupyter`
1. Start changing the code as described in the notebook.
   The instructions will ask you to create two new files in the same folder `source.py` and `model_specific_helpers_2.py`
   The example contains instructions to outsource the code into different those files and streamline the notebook.
1. Although the example is very long since it covers ALL the steps that were necessary for the visit model, YOUR version
   should be very short in the end. (Like instpectModel.py) 

### General Remarks and Additional information

## Keep your repository updated! ## (Ohterwise you will base your changes on outdated code run into conflicts and make it very hard for us to incorporate your changes (and for you to incorporate ours...) 

1. If you cloned the main repository ( https://github.com/MPIBGC-TEE/bgc_md2 ) it's a good practice to 
   say `git pull` from time to time and definitely before you start working and change files.
  
1. If you use your own fork you should also update from the upstream repository frequently.
   The first time you do this you have to tell your git which is the upstream repo **just once**:
   ```
   cd into/cloned/fork-repo
   git remote add upstream https://github.com/MPIBGC-TEE/bgc_md2
   git fetch upstream
   ```
   every time you want to pull the latest changes into your master you say:
   ```
   git pull upstream master
   ```
   
## Don't push directly to the master of the main repo ##
- If you have your own fork you are already save.
- In case you cloned the main repository the simplest strategy is to work with a branch.
  https://git-scm.com/book/en/v2/Git-Branching-Branches-in-a-Nutshell
  
  - create a branch of your own (the `-b` is only necessary the first time)
    ```
    git checkout -b   {YourUniqueBranchName} 
    ```
  - **or** use the `test` branch 
    ```
    git checkout test
    ```  
  and after committing your changes 
  push to it's remote.
  This will trigger a run of bgc_md2's testsuite on github, regardless of the branch you used.
  You can see the results under the `Actions` tab on the github page https://github.com/MPIBGC-TEE/bgc_md2/actions.  
  **Only if everything works** we merge the commits of your branch into the master and push again... 
  ```
  git checkout master 
  git merge {YourUniqueBranchName}
  git push
  ```
  Thus the master branch of the main repo stays reliable for everybody and *failed/incomplete experiments* 
  (including unintended ones)  do not affect your colleagues.

 
## Run the tests 
To do the actual scientific work of comparing different models we have to be able to call some of your model specific functions with some specific arguments. E.g. we will use your optimised parameters in a tracability plot.
To make sure that we can do so for your model cd into your folder and run 
```
python ../run_my_tests.py
```
This will execute the test suite ```TestSymbolic.py``` for your model. We run it for all models regularly
but it's best if you notice problems with your own model first.
If all the tests pass, we can use your model for model inter comparisons. 
As we extend the analysis more tests will be added which means that you might have to add or change a function to `model_specific_helpers_2.py` or a parameter to `test_helpers.py` 
There will however usually be already an example model available that passes the tests.
The commit messages (you can display the last 4 with `git log -4`) give you a hint about the purpose of the change and `git log -p -4` shows you the code changes of the last 4 commits.
You can find more tricks here:https://git-scm.com/book/en/v2/Git-Basics-Viewing-the-Commit-History

Think of the test suite as an insurance against unwelcome suprises and a way to make it easier for others to look at your results or help you with a specific problem. 
Run it often, at least once before you check in!
If you want to run a single test (e.g. because the test `test_symobolic_description` fails) you can run:
```
python ../run_my_tests.py test_symobolic_description
```
Make sure that you pull changes often (and allways before you start editing) because we now really work on several models at the same time, so many small things will change often.





