Structure
=========
 * Schema
   The package represents domain knowledge by two sets:

   * A set of variables called `MVars` and 
   * A set of functions called `Computers`, operating on sets of instaces of these types.

   Sets of `MVars` form the nodes and the computers the edges of a directed graph G.
   We call G the computability graph because a desired `Mvar` (target) can be computed from 
   a set of `Mvars` (src) if there is a path from src to target.

   `MVars` are instanciated and combined into `MVarSets` by user provided python code, typically placed in a file `models/NameOfTheModel/source.py`.
   This is however just a convention (with some support by helper functions). It is also possible to place the code elsewhere.
   The `models` submodule of `bgc_md2` can thus be seen as a set of examples
   The (minimal) purpose of the model file is usually that an instance of 
   either a CompartmentalModel or a ModelRun can be created.
   This can be done
   * directly, by calling one of the constructors (of a CompartmentalModels) or 

   * indirectly by some definitions of variables with special names from which the ingredients for the constructor calls can be infered.

   The set of computers represented by the `bgc_md/resolve/MvarsAndComputers.py` 


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
     
     *  Viewed from a high level of abstraction (symbolic) models and (climatic driver) data seem to be orthogonal to each other. 
        but they are also connected by a mapping between the dataset and the models realized by model specific data preparation code.
        A model with more pools and connections we will
        need more data and will connect the data to different pools or fluxes. 
        The data representation code  will reference one model but will be used for different (similar datasets).
