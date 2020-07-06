## Workflows To develop a useinterface for the framework it is essential to
know the applications.

* For administrators (group leader) it is important to produce an overview of
  the work that has been done with the framework. 
  This has been achieved up to now by a static website created from
  the yaml files describing the models and report templates (recursively
  defined snippets) for output creation.
  * With the new frameork it should still be possible to get an overview quikly
    (Overview Notebook)
  * It would be preferable if the Notebook would reflect changes in the models
    immediately (dynamically) instead of having to wait for the `compilation` of a
    report.
  * The Overview page of the static website linked to model specific reports.
    These very much depend on the information available and can in some cases
    involve very expensive computations. As a result the creation of the whole
    static website could take a long time.  It would be much better if the new
    framework would do these computations on demand, but cache them.  Since
    jupyter notebooks contain their output it should be possible to keep a
    collection of available reports on github (or any other site offering to
    host *.ipynb files 

* For model investigators (including first and foremost the researchers of the group)
  it is important to be able to define new models wiht ease to get the results quickly.
  There are two main objectives.
  * use (and create) the common infrastructure. For a large part this is achieved by the 
    use of CompartmentalSystems and LAPM 
  * To compare models (mainly Veronika) the consistent representation of modelproperties is important 
    (The main reason to have this framework)
    Part of this work is the translation of the yaml files. The process could be assisted by showing computable
    properties.

  * for new users it should be possible 
    * start the development of new models quikly by copying existent prototypes
    * explore the possibilities for analysis (This is where the computational graph should be of assistence)


   
