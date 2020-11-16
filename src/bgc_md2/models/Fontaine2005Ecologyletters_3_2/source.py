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
        'C_s': 'carbon stock in soil organic matter' # "\\text{quantitiy of carbon}"
        ,'C_f': 'carbon stock in fresh organic matter' # "\\text{quantitiy of carbon}"
        ,'C_ds': 'carbon stock in soil organic matter decomposers' # "\\text{quantitiy of carbon}"
        ,'N': '"mineral nitrogen pool ($N:C$ ratio always constant)"' # "\\text{quantitiy of nitrogen}"
        ,'A': 'decomposers consumption rate of SOM' # "\\text{time}^{-1}"
        ,'r': '"fraction of decomposer biomass released as CO$_2$"' # "\\text{time}^{-1}"
        ,'s': '"decomposers production rate of soil organic matter"' # "\\text{time}^{-1}"
        ,'k': '"rate of fresh organic matter decomposition under substrate limitation ($N$ excess)"' # "\\text{time}^{-1}"
        ,'alpha': '"$N:C$ ratio in soil organic matter and in decomposers"'
        ,'beta': '"$N:C$ ratio in fresh organic matter"'
        ,'i': '"rate of mineral $N$ diffusion in soil"' # "\\text{time}^{-1}"
        ,'Phi_l': 'fresh organic matter carbon flux' # "(\\text{quantity of carbon})(\\text{time}))^{-1}"
        ,'Phi_i': 'nitrogen that flows into the ecosystem' # "(\\text{quantity of nitrogen})(\\text{time}))^{-1}"
        ,'Phi_o': 'nitrogen that flows out of the ecosystem' # "(\\text{quantity of nitrogen})(\\text{time}))^{-1}"
        ,'Phi_up': 'nitrogen flux associated with the nitrogen uptake by the plant cover' # "(\\text{quantity of nitrogen})(\\text{time}))^{-1}"
}

for name in sym_dict.keys():
    var(name)
t = TimeSymbol("t") # unit: ""
x = StateVariableTuple((C_s, C_f, C_ds, N))
u = InputTuple((0, Phi_l, 0, Phi_i-Phi_o-Phi_up))
B = CompartmentalMatrix([[0,              0,                        s-A,               0],
                         [0,              0,      -alpha*r/(alpha-beta), -i/(alpha-beta)],
                         [0,              0, A-s-r+alpha*r/(alpha-beta),  i/(alpha-beta)],
                         [0,              0,                          0,              -i]])

mvs = MVarSet({
    BibInfo(# Bibliographical Information
        name="FB2005 - The C-N model of SOM dynamics with a single type of decomposer (model 3)",
        longName="The C-N model of SOM dynamics, two decomposer types", 
        version="2 - N limitation",
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
})
