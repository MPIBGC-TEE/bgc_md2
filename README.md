
## Installation
```bash 
source installation_developer_conda.sh
```
## Structure
This prototype represents models as subfolders of folder 
```bash
models
```
Every model folder contains at least one file 
```bash
source.py
```
That can contain import other module that can live in the same directory or be accessible globally. 
The (minimal) purpose of the model file is usually that an instance of either a CompartmentalModel or a ModelRung can be created.
This can be done
* directly, by calling one of the constructors (of a CompartmentalModels) or 
* indirectly by some definitions of variables with special names from which the ingredients for the constructor calls can be infered.

The package contains (growing) knowledge about these inferences. 
This knowledge is represented in two sets.
A set of variables called `MVars` and a set of functions called `Computers`.
The sets form the nodes and edges of a directed graph, which can be used to 
determine if a desired `Mvar` can be computed from the `Mvars` defined (or instanciated) in the model descritpion in `models/NameOfTheModel/source.py` by the users and the (growing) set of computers represented by the `bgc_md/resolve/MvarsAndComputers.py` 

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


## Ensurance of reaching the objectives
The development should be test driven. The required functionality is represented by the 'reports' we want to create, which intuitively translates into a set of jupyter notebooks we want to keep functional regardless of changes in the
framework or the model `source.py` files.
The set of functional reports may contain notebooks of both kinds:
* model specific analyses (plural)
* comparisons 
   * of models
   * of model runs (for numeric cases where the models may not be inferable from data)
It is evident that it is undesirable to maintain duplicated similar notebooks for different models.
Either the (report=notebook) creation has to be automated or the availability of common results for sets of models has to be tested. This shows that the ability to compare models is a byproduct of deduplication.
From a testing perspective a set of notebooks is unwieldy though. It has the role of integration tests, which are much more complex and maintainance intensive that unit test. 
It seems more reasonable to create a test matrix (as the tensor product set of the set of request and a set of models ).  


## various notes on implementation

* The 'Computers' and 'MVars' represent a set of types and strictly typed
  functions (including the return values) This could be (re)implemented with
  the new python type annotations.  An advantage is that we can express our
  idea in a well defined and well documented way and avoid extra effort for the
  user..  
* The computibility graph is expensive to create and only changes if new
  `Computers` and `MVars` are created.  It should be cached, which encurages
  the use of immutable data structures. (since we can use functools )



## to do:
* With a growing set of computers it is possible that some of the MVars
  explicitly defined by a user in a specific model could be computed from a
  subset.  This introduces an indirect form of duplication, which should be
  detected automatically and flagged (so that the model descriptions can be
  pruned to a minimal set)  Even more indirectly some results may be computable
  in different ways from the intanciated MVars (some MVars have several edges
  pointing to them) Contradictions should be checked for.

* Library vs. Framework CompartmentalSystems and LAPM can be considered
  libraries that can be called in different parallelisation scenarios and with
  different data.  
  bgc-md has both library and framework characteristics.  
   * 	The (types) of MVars and computers connecting them represents a domain knowledge
    	in form of a library.  
   * 	The collection of different models enforces some rules on the way a
     	single model has to be formulated in order to be comparable to other
	models. 
	In this sense it comprises a framework into which different models can be
	plugged.  
  This raises practical questions.
  (Symbolic) Models and data can (preferably) be seen as orthogonal to each other, 
  but are also connected. 
  (for a model with more pools and connections we will
  need more data...) 
  Since the aim is to compare models across different datasets the model specific part of 
  data representation should be minimized.
  How much (links to the) data should be part of the database?

 * Translation of yaml files to the new format.
   * examples: 
     * ``` bgc_md2/bgc_md2/models/Williams2005GCB/source.py ``` 
     * ``` bgc_md2/bgc_md2/models/Potter1993GlobalBiogeochemicalCycles/source.py```
   * In the yaml file find all variables that do not have an "expr" (are not defined as expressions with respect to other variables)
     create ```Symbol``` Instances for them.
     (This can be done at the same time as providing a bibliographic variable dictionary which is only used for presentation)
     If no parametrisations are provided the dimensions can be kept as comments.
     If quantitative model runs are desired the units are actually used for the values of those variables so they 
     belong to the parameter sets.
     The units like ```kilogramm```,```meter``` have to be imported from sympy.physics.units 
   * Move the remaining variables which are defined by expression into a section after this and copy the decription as a comment.
     The dimension  does not have to be set for these because it can be inferred from the dimensions of the variables
     constituting the expression.
   * In the model file look out for common variables (present in more than one model) that are already defined as
     Classes in ```bgc_md2/resolve/mvars.py```
     * Import the appropriate constructors in the ```source.py``` file and use them directly for the variables that go into
       the ```specialVars``` set. These are the things that we can compare between models that define them 
       (or other specialVars from which they can be derived)
     * If you want a variable be comparable across models and you do not find it in ```bgc_md2/resolve/mvars.py``` yet, then
       note it down in this document in the following bullet list, like the first point and give an example of a succesfull comparison.
	 *  mvars to be implemented:
	  * Implement a notebook wich compares vegetation cycling matrices vegetation pools across models.
  	    steps:
	    * This could internally be done via the dyadic representation (but should also be possible if the model is given in matrix form) 
  	    * Actually the stateVariableTuple should be given externally (independent from the one in the mode) so that
	      the same order can be used to be able to compare the blockstructure.
	    * find models described by the old yaml files that have a similar or identical structure for the vegetation sub matrix if similar coordinates are used.
	      probably a unified nomenclature for Leaf Root and Labile pools has to be used to make them comparable.
	    * How much of this can be (at first) implemented outside tme mvars and computers modules just in the notebook or python file where the comparison is made?
	    * build a test.

* the dependency list for conda is probably to full (just copied from the widget tutorial)
  After the UI is more or less stable the relavant packages have to be isolated and should also be moved into setup.py if no version pinning is required)

