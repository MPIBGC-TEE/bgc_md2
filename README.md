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
  idea in a well defined and well documented and avoid extra effort for the
  user..  The challenge would be to perform the computebility inference
  (including the assembly of defined Mvars=types) at runtime (from inside a
  jupyter notebook). At the moment it is unclear which tools or libraries we
  can use for this purpose.  
* The computibility graph is expensive to create and only changes if new
  `Computers` and `MVars` are created.  It should be cached, which encurages
  the use of immutable data structures. (since we can use functools )



## current problems 
* At the moment the tests can not be run in one suite (we have to find the
  reason)a, since the source of the problem most likely lies in the
  (prototypically implemented) way the package handles more than one model and
  could affect model comparisons.

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
  (Sympolic) Models and data can (preferably) be seen as orthogonal to each other, 
  but are also connected. 
  ( for a model with more pools and connections we will
  need more data...) 
  Since the aim is to compare models across different datasets the model specific part of 
  data representation should be minimized.
  How much (links to the) data should be part of the database?


