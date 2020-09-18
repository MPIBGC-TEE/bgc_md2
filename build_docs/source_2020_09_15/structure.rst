Structure
=========
 * Schema
   The package represents domain knowledge by two sets:
   * A set of variables called `MVars` and 
   * A set of functions called `Computers`, operating on sets of instaces of these types.

   The instance sets form the nodes and the computers the edges of a directed graph G.
   We call G the computability graph because a desired `Mvar` (target) can be computed from 
   a set of `Mvars` (src) if there is a path from src to target.
   We can also use

   defined (or instanciated) in the user provided model descritpion in `models/NameOfTheModel/source.py` by the users and the (growing) set of computers represented by the `bgc_md/resolve/MvarsAndComputers.py` 
Every model folder contains at least one file 
```bash
source.py
```
That can contain import other module that can live in the same directory or be accessible globally. 
The (minimal) purpose of the model file is usually that an instance of either a CompartmentalModel or a ModelRun can be created.
This can be done
* directly, by calling one of the constructors (of a CompartmentalModels) or 
* indirectly by some definitions of variables with special names from which the ingredients for the constructor calls can be infered.
   
   The package does not use standard database technology like SQL but stores available data connected to a model as a set of 
   variables of predefined types (from now on called MVars).
   In terms of a database one can think of every Mvar type as a table  
   with two columns model_id 
 * Models
 * Notebooks
   * For single models
   * For model comparisons

 * Notes on implementation
   * Library vs. Framework 
     CompartmentalSystems and LAPM can be considered
     libraries that can be called in different parallelisation scenarios and with
     different data.  
     bgc-md has both library and framework characteristics.  
     * 	The (types) of MVars and computers (functions) connecting them represents a domain knowledge
       	in form of a library.  
     * 	The collection of different models enforces some rules on the way a
        single model has to be formulated in order to be comparable to other
        models. 
        In this sense it comprises a framework into which different models can be
        plugged.  
     
     *  Viewed from a high level of abstraction models (symbolic) models and (climatic driver) data seem to be orthogonal to each other. 
        but they are also connected by a mapping between the dataset and the models realized by model specific data preparation code.
        A model with more pools and connections we will
        need more data and will connect the data to different pools or fluxes. 
        The data representation code  will reference one model but will be used for different (similar datasets).
