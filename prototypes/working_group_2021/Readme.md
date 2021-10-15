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
         git clone https://github.com/MPIBGC-TEE/bgc_md2.git
         ``` 
         in your terminal 
         (in the directory where you want to checkout the code)
      
## How to use the code
### Structure of the directory
The files are organized in two catagories:
 * model specific
 * general
For every model there is a specific sub folder, containing model the specific files, while the general functions reside one level above.
The first model (and at the time of writing the only) represented in this way is the modified version of  Yuanyuan's cable example.
The folder name is  `yy_cable`.
### Test the example
Before you start working on your own version of it you should make sure that you can run the example.
1. Depending on your python environment you will have to install some packages: 
  * json
  * numpy
  * pandas
  * netCDF4
  * tqdm
1. `cd yy_cable`Create a small config file named `config.json` to tell the code where to look for the data: 
In my case the files reside in `/home/data/yuanyuan` and the content of `config.json` is 
```json
{"dataPath": "/home/data/yuanyuan"}
```
1. After the successful installation you can change into `yy_cable` and should be able to execute the mcmc by typing`python main.py`

### Adapt the code to your own needs
1. copy and rename the folder 
1. start adapting the code to your own model
   There are several ways how to do this depending on your workflow.
   I what follows I describe a how I went about adapting Alisons code.
   I use automated tests to get one part after the other to work and conserver the progress.

  To do this it step by step
	1. run the TestSuite in the newly created folder by 
           ```bash
           python -m unittest TestModel.py 
	   '''
           To make sure that you start with an example that works.

	1. Open the files  `TestModel.py' and `model_specific_helpers.py'
	   and change the names and values of the `UnEstimatedParameters' and `EstimatedParameters'
	   in the definition of `model_specific_helpers.py' 


### General Remarks and Additional information


