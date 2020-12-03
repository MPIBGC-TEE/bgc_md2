import numpy as np
from sympy import var
from frozendict import frozendict
from bgc_md2.resolve.mvars import (
    CompartmentalMatrix,
    InputTuple,
    TimeSymbol,
    StateVariableTuple,
)
from ..BibInfo import BibInfo 
#from bgc_md2.resolve.MVarSet import MVarSet
from bgc_md2.helper import MVarSet

sym_dict = {
        'A': 'labile pool' # "MgC*ha^{-1}" 
        ,'B': 'stable pool' # "MgC*ha^{-1}" 
        ,'alpha': 'annual decomposition rate of labile pool' # "yr^{-1}"
        ,'beta': 'annual decomposition rate of stable pool' # "yr^{-1}"
        ,'m': 'annual organic matter input' # "MgC yr^{-1}"
        ,'K': 'isohumic coefficient'
}

for name in sym_dict.keys():
    var(name)
t = TimeSymbol("t") # unit: "year"
C = StateVariableTuple((A,B))
I = InputTuple((m,0))
A_GeM = CompartmentalMatrix(
[[-alpha,     0 ],
 [alpha*K, -beta]])

mvs = MVarSet({
    BibInfo(# Bibliographical Information
        name="",
        longName="Henin and Dupois", 
        version="1",
        entryAuthor="Holger Metzler",
        entryAuthorOrcid="0000-0002-8239-1601",
        entryCreationDate="09/03/2016",
        doi="",
#        bibtex="@inproceedings{Henin1945Annalesagronomiques,
#                     author = {H\\'{e}nin, S and Dupuis, M},
#                     booktitle = {Annales agronomiques},
#                     pages = {17--29},
#                     title = {{Essai de bilan de la mati\\`{e}re organique du sol}},
#                     volume = {15},
#                     year = {1945}
#                    }",
        sym_dict=sym_dict
    ),
    A_GeM,  # the overall compartmental matrix, decomposition operator
    I,  # the overall input
    t,  # time for the complete system
    C,  # state vector of the complete system
})


