from sympy import var
from frozendict import frozendict
from bgc_md2.resolve.mvars import (
    CompartmentalMatrix,
    InputTuple,
    TimeSymbol,
    StateVariableTuple,
)
from ..BibInfo import BibInfo 
from ComputabilityGraphs.CMTVS import CMTVS
from bgc_md2.helper import bgc_md2_computers

sym_dict = {
        'C_s': 'carbon stock in soil organic matter' # "\\text{quantitiy of carbon}"
        ,'C_ds': 'carbon stock in soil organic matter decomposers' # "\\text{quantitiy of carbon}"
        ,'A': 'decomposers consumption rate of SOM' # "\\text{time}^{-1}"
        ,'r': '"fraction of decomposer biomass released as CO$_2$"' # "\\text{time}^{-1}"
        ,'Phi_l': 'fresh organic matter carbon flux' # "(\\text{quantity of carbon})(\\text{time}))^{-1}"
}

for name in sym_dict.keys():
    var(name)
t = TimeSymbol("t") # unit: ""
x = StateVariableTuple((C_s, C_ds))
u = InputTuple((0, Phi_l))
B = CompartmentalMatrix([[0,  -A],
                         [0, A-r]])

mvs = CMTVS(
    {
        BibInfo(# Bibliographical Information
            name="FB2005 - The model of SOM decomposition (model 1)",
            longName="The C-N model of SOM dynamics, two decomposer types", 
            version="1",
            entryAuthor="Holger Metzler",
            entryAuthorOrcid="0000-0002-8239-1601",
            entryCreationDate="22/03/2016",
            doi="10.1111/j.1461-0248.2005.00813.x",
            sym_dict=sym_dict
        ),
        B,  # the overall compartmental matrix
        u,  # the overall input
        t,  # time for the complete system
        x,  # state vector of the complete system
    },
    bgc_md2_computers()
)
