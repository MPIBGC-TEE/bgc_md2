from sympy import var, ImmutableMatrix, diag
from bgc_md2.resolve.mvars import (
    CompartmentalMatrix,
    InputTuple,
    TimeSymbol,
    StateVariableTuple,
    VegetationCarbonInputScalar,
    VegetationCarbonInputPartitioningTuple,
    VegetationCarbonStateVariableTuple,
)
from ..BibInfo import BibInfo 
from ComputabilityGraphs.CMTVS import CMTVS
from bgc_md2.helper import bgc_md2_computers

sym_dict = {
        'C_il':      'Carbon in leaves of plant functional type (PFT) $i$',
        'C_is':      'Carbon in transport tissue (mainly stems) of PFT$_i$',
        'C_ir':      'Carbon in fine roots of PFT$_i$',
        'NPP_i': 'Net Primary Production for PFT$_i$',
        'a_il':  'Fraction of annual NPP allocated to leaves for PFT$_i$',
        'a_is':  'Fraction of annual NPP allocated to stem for PFT$_i$',
        'a_ir':  'Fraction of annual NPP allocated to roots for PFT$_i$',
        'tau_il': 'Residence time of carbon in leaves for PFT$_i$  ',
        'tau_is': 'Residence time of carbon in stem for PFT$_i$  ',
        'tau_ir': 'Residence time of carbon in roots for PFT$_i$ ',
        'S': 'Percent sand in soil',
}

for name in sym_dict.keys():
    var(name)

a_il = -0.0025 * S + 0.44 
a_ir = 0.0039 * S + 0.137
a_is = 1 - a_il - a_ir

x = StateVariableTuple((C_il, C_is, C_ir))
u = NPP_i
b = (a_il, a_is, a_ir)
Input = InputTuple(u*ImmutableMatrix(b))
A = CompartmentalMatrix(
    diag(-1/tau_il, -1/tau_is, -1/tau_ir)
)
t = TimeSymbol("t")

# The following parameter set corresponds to a previous version of the model with fixed coefficients
#model_run_data:
#    parameter_sets:
#        - "Tropical evergreen trees":
#            values: {a_il: 0.25, a_is: 0.5, a_ir: 0.25}
#            doi: 10.5194/bg-10-2255-2013

mvs = CMTVS(
    {
        BibInfo(# Bibliographical Information
            name="IBIS",
            longName="Integrated Biosphere Simulator", 
            version="2.6",
            entryAuthor="Verónika Ceballos-Núñez",
            entryAuthorOrcid="0000-0002-0046-1160",
            entryCreationDate="22/3/2016",
            doi="10.5194/bg-10-2255-2013",
            #further_references=BibInfo(doi="10.5194/bg-10-2255-2013"),
            sym_dict=sym_dict
        ),
        A,  # the overall compartmental matrix
        Input,  # the overall input
        t,  # time for the complete system
        x,  # state vector of the complete system
        VegetationCarbonInputScalar(u),
        # vegetation carbon partitioning.
        VegetationCarbonInputPartitioningTuple(b),
        VegetationCarbonStateVariableTuple((C_il, C_is, C_ir)),
    },
    bgc_md2_computers()
)
