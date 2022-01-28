# Unified data assimilation resources for the working group
## How to get the code and contribute:
This folder is part of a bigger repository.
To use the code You have two options.
1. Fork the repository:
    * Pro:
      Forking (creating your own branch of the repo) has the advantage that it is harder to mess up anybody else's work,
      even if you are not too familiar with git, since you can not directly push to the master.
    * Con:
      To make your changes available to the group you have to make a pull request and I will pull your code.
    * How to
      1. read the general github instructions here: https://docs.github.com/en/get-started/quickstart/fork-a-repo 
      1. press the Fork button at the upper right cornder of your github window.
      1. Clone you own fork
2. Clone the repository 
   * Pro:
      I will grant you the priveliges to push which allows you to make changes available quicker to everybody.
    * Con:
      If you don't know what a merge conflict is or how to resolve it without messing up somebody else's work then we are all better off 
      with you requesting us to pull your changes...
    * How to: 
      1. Go to the main site of the repo: https://github.com/MPIBGC-TEE/bgc_md2
      2. Press the green Code button and copy the url
      3. type 
         ```bash
         git clone https://github.com/MPIBGC-TEE/bgc_md2.git --recurse submodules
         ``` 
         in your terminal 
         (in the directory where you want to checkout the code)
      
## How to use the code
### Structure of the directory
The files are organized in two catagories:
 * model specific
 * general
For every model there is a specific sub folder, containing the model specific files, while the general functions reside one level above.
The first model  represented in this way is the modified version of  Yuanyuan's cable example.
The folder name is  `yy_cable`.
### Test the example
Before you start working on your own version of it you should make sure that you can run the example.
1. Install the `bgc_md2` package as described in the instructions in the gitHub https://github.com/MPIBGC-TEE/bgc_md2#installation
   Remark: If you cloned the repo change the url to your fork, otherwise you will not see your latest contributions.

1. `cd yy_cable`Create a small config file named `config.json` to tell the code where to download and look for the data: 
In my case the files reside in `/home/data/yuanyuan` and the content of `config.json` is 
```json
{"dataPath": "/home/data/yuanyuan"}
```
1. Run the little script to download the data.
1. After the successful installation you can change into `yy_cable` and should be able to execute the mcmc by typing`python main.py`

### Adapt the code to your own needs
1. Create a folder for your model. We call it `{your_model}` from now on.
1. Start copying and changing the code and the automated tests for your own model
   To do this step by step
   	1. Create a symbolic description of the model.
	   Check the paper that describes your model and find out if the equations are already in matrix form
	   or if they are given as fluxes (into, out of , or between pools). 
	   This decides which example you should start with.
	   An example folder for a matrix description is `yy_cable`
	   The folder `kv_visit2` contains an example starting from fluxes.
	   The package will create the matrix from this input for you.
	   We will call the folder you decide to start with `{example}` from now on.  

   		1. copy the following files from `{example}` into your new folder:
		   ```bash
		   run_tests_serial.py, TestSymbolic_1.py` into your new folder
		1. run the TestSuite in the newly created folder by 
        	   ```bash
        	   python -m unittest TestModel.py 
		   '''
        	   To make sure that you start with an example that works.

	1. Open the files  `TestModel.py' and `model_specific_helpers.py'
	   and change the names and values of the `UnEstimatedParameters' and `EstimatedParameters'
	   in the definition of `model_specific_helpers.py' 


### General Remarks and Additional information


