import numpy as np
from sympy import var, Symbol, ImmutableMatrix, exp
from frozendict import frozendict
from bgc_md2.resolve.mvars import (
    CompartmentalMatrix,
    InputTuple,
    TimeSymbol,
    StateVariableTuple,
    NumericParameterization,
    NumericStartValueDict,
    NumericSimulationTimes,
)
from ..BibInfo import BibInfo
from ComputabilityGraphs.CMTVS import CMTVS
from bgc_md2.helper import bgc_md2_computers

sym_dict = {
        'X': 'microbial biomass pool',  # "\\mu gC cm^{-3}"
        'S': 'substrate pool',  # "\\mu gC cm^{-3}"
        'mu_max': 'maximal relative growth rate of bacteria',  # "hr^{-1}"
        'D_max': 'maximal relative death rate of bacteria',  # "hr^{-1}"
        'K_s': 'substrate constant for growth',  # "\\mu gC cm^{-3}"
        'K_d': 'substrate constant for death of bacteria',  # "\\mu gC cm^{-3}"
        'K_r': 'fraction of dead biomass recycling to substrate',
        'theta': 'soil water content',
        # "ml\\text{ solution }cm^{-3}\\text{ soil}"
        #
        'Y': 'yield coefficient for bacteria',
        'ExuM': 'maximal exudation rate',  # "\\mu gC hr^{-1}cm^{-3}"
        'ExuT': """time constant for exudation, 
                responsible for duration of exudation""",  # "hr^{-1}"
        'BGF': 'constant background flux of substrate',  # "\\mu g C cm^{-3}hr^{-1}"
        'mu_S': 'relative growth rate of bacteria (dependent on substrate concentration)',
        'Exu': 'exudation rate (dependent on time)',
}
for name in sym_dict.keys():
    var(name)

# t is not an instance of  symbol but of TimeSymbol, Therefore it can not
# created along with the other symbols but hast to be created before it appears
# in  any expressions.
t = TimeSymbol('t')  # in the original pub unit: "hour"

mu_S = mu_max * S/(K_s*theta+S)  # "hr^{-1}"
Exu = ExuM * exp(-ExuT*t)  # "hr^{-1}"
x = StateVariableTuple((X, S))
u = InputTuple((0, BGF + Exu))
T = ImmutableMatrix([[ -1,  Y],
                     [K_r, -1]])
N = ImmutableMatrix([[D_max*K_d/(K_d+S/theta),                        0],
                     [                      0, X/Y*mu_max/(K_s*theta+S)]])
B = CompartmentalMatrix(T * N)
# Original values from linked model (no nitrogen cycle considered in this model here):
np1 = NumericParameterization(
    par_dict={
D_max: 0.26, Y: 0.44, K_s: 3.0, theta: 0.23, ExuT: 0.8, BGF: 0.15, K_d: 14.5, K_r: 0.4, ExuM: 8.0, mu_max: 0.063},
    func_dict=frozendict({})
)

# Standard version of BACWAVE, optimized to simulate bacterial biomass along wheat roots:
nsv1 = NumericStartValueDict({X: 0.5, S: 1.5}) # "Low"
nsv2 = NumericStartValueDict({X:1.0, S: 2.5}) # "Medium"
nsv3 = NumericStartValueDict({X: 1.5, S: 4.0}) # "High"
ntimes = NumericSimulationTimes(np.arange(0,2000,0.1))

mvs = CMTVS(
    {
        BibInfo(# Bibliographical Information
            name="BACWAVE",
            longName="", 
            version="",
            entryAuthor="Holger Metzler",
            entryAuthorOrcid="0000-0002-8239-1601",
            entryCreationDate="15/03/2016",
            doi="10.2307/4251775",
            sym_dict=sym_dict
        ),
        B,  # the overall compartmental matrix
        u,  # the overall input
        t,  # time for the complete system
        x,  # state vector of the complete system
        np1,
        nsv1,
        ntimes
    },
    bgc_md2_computers()
)
