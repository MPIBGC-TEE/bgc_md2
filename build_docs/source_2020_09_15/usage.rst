Usecases / Objectives
=====================
* Exploration of existing models:
  The package is supposed to assist in the creation of 'reports' in the form of jupyter notebooks.
  The notebooks will be of the following types:
  * Investigations of a single model (or modelrun).

  * Comparisons between models (modelruns).

* Assistence in the formulation of new models in a step by step process.
  To this end the framework implements some powerfull interactive tools.
  In this scenario the role of the `bgc_md` package is to guide the user (=
  author of a particular notebook concerned with a particular model, and
  simultaniously author of the `source.py` of that model) by using the
  computability graph (as represented by `bgc_md/resolve/MvarsAndComputers.py`)
  to either 
  * show which addidional results can be computed, given the
    information already present in the models `source.py` or 
  * show which additional information has to be provided in the models `source.py` 
    to be able to obtain a desired result.

