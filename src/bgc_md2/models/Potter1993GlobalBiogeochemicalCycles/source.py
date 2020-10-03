from sympy import var, symbols, Symbol, ImmutableMatrix, diag, exp, Rational
from frozendict import frozendict
from bgc_md2.resolve.mvars import (
    CompartmentalMatrix,
    InputTuple,
    TimeSymbol,
    StateVariableTuple,
    VegetationCarbonInputScalar,
    VegetationCarbonInputPartitioningTuple,
    VegetationCarbonStateVariableTuple,
    NumericParameterization,
    NumericStartValueDict,
    NumericParameterizedSmoothReservoirModel,
    NumericSimulationTimes,
)
from ..BibInfo import BibInfo 
from bgc_md2.helper import MVarSet

# We can retain make the original description available to the framework as 
# part of the bibliographical information in a dictionary which could also be used to define
# the symbols (to avoid duplication) as demonstrated here.
sym_dict={
        'C_f':      'Carbon in foliage',
        'C_r':      'Carbon in roots',
        'C_w':      'Carbon in woody tissue',
        'SOL':      'Total solar radiation (SOL(x,t))',
        'epsilon':  """PAR use efficiency ($\\epsilon(x,t)$). Function that
                    depends on effects of temperature and water stress""",
        'FPAR':     """ Fraction of incoming PAR intercerpted 
                    by green vegetation (FPAR(x,t))""",
        'IPAR':     """Intercepted photosynthetically active
                    radiation(IPAR(x,t)).  The factor of 0.5 accounts for the
                    fact that approx. half of SOL is in PAR waveband (0.4-0.7
                    $\\mu$m)""",
        'alpha_f':  """Proportional allocation constant of available carbon
                    allocated to foliage""",
        'alpha_r':  """Proportional allocation constant of available carbon
                    allocated to roots""",
        'alpha_w':  """Proportional allocation constant of available carbon
                    allocated to wood""",
        'tau_f':    'Inverse Residence time of carbon in foliage',
        'tau_r':    'Inverse Residence time of carbon in roots',
        'tau_w':    'Inverse Residence time of carbon in wood'
}
for name in sym_dict.keys():
    var(name)

IPAR = SOL*FPAR*0.5
NPP = IPAR * epsilon
x = StateVariableTuple((C_f, C_r, C_w ))
u = NPP
b = (alpha_f, alpha_r, alpha_w)

Input = InputTuple(u * ImmutableMatrix(b))
A = CompartmentalMatrix(
    diag(-tau_f, -tau_r, -tau_w)
)
t = TimeSymbol("t")

# model_run_data is incomplete :
# Even parameter_sets cannot fully parameterize the symbolic model: mm made up
# a purely fictional parameters to show how a purely numeric could be build.
# (The yaml file lacks consistent unit information).
np1 = NumericParameterization(
    par_dict={
        alpha_f:    Rational(1, 3),
        alpha_r:    Rational(1, 3),
        alpha_w:    Rational(1, 3),
        tau_f:      1.5,
        tau_r:      3,
        tau_w:      50,
        SOL:        1, # completely faked by mm
        epsilon:    1.1,
        FPAR:       1, # completely faked by mm
    },
    func_dict=frozendict({})
)
# Original dataset of the publication":
#            values: {alpha_f: 'Rational(1,3)',alpha_r: 'Rational(1,3)',alpha_w: 'Rational(1,3)'}
# # epsilon varies from 1.1 - 1.4 gC*MJ^{-1}PAR in crop ecosystems.
#            doi: 10.1029/93GB02725
#        - "Tundra":
#            values: {alpha_f: 0.25, alpha_r: 0.25, alpha_w: 0.5, tau_f: 1.5, tau_r: 3, tau_w: 50}
#            doi: 10.2307/1313568
#        - "High-latitude forest":
#            values: {alpha_f: 0.30, alpha_r: 0.25, alpha_w: 0.45, tau_f: 1, tau_r: 3, tau_w: 50}
#            doi: 10.2307/1313568
#        - "Boreal coniferous forest":
#            values: {alpha_f: 0.25, alpha_r: 0.25, alpha_w: 0.5, tau_f: 2.5, tau_r: 3, tau_w: 50}
#            doi: 10.2307/1313568
#        - "Temperate grassland":
#            values: {alpha_f: 0.45, alpha_r: 0.55, tau_f: 1.5, tau_r: 5}
#            doi: 10.2307/1313568
#        - "Mixed coniferous forest":
#            values: {alpha_f: 0.25, alpha_r: 0.25, alpha_w: 0.5, tau_f: 1.5, tau_r: 3, tau_w: 40}
#            doi: 10.2307/1313568
#        - "Temperate deciduous forest":
#            values: {alpha_f: 0.30, alpha_r: 0.25, alpha_w: 0.45, tau_f: 1, tau_r: 3, tau_w: 40}
#            doi: 10.2307/1313568
#        - "Desert and bare ground":
#            values: {alpha_f: 0.25, alpha_r: 0.25, alpha_w: 0.5, tau_f: 1.5, tau_r: 3, tau_w: 50}
#            doi: 10.2307/1313568
#        - "Semi-arid shrubland":
#            values: {alpha_f: 0.25, alpha_r: 0.25, alpha_w: 0.5, tau_f: 1.5, tau_r: 3, tau_w: 50}
#            doi: 10.2307/1313568
#        - "Savanna and woody grassland":
#            values: {alpha_f: 0.30, alpha_r: 0.25, alpha_w: 0.45, tau_f: 1, tau_r: 5, tau_w: 25}
#            doi: 10.2307/1313568
#        - "Tropical evergreen rain forest":
#            values: {alpha_f: 0.25, alpha_r: 0.25, alpha_w: 0.5, tau_f: 1.5, tau_r: 2, tau_w: 25}
#            doi: 10.2307/1313568
#specialVars = {
mvs=MVarSet({
    BibInfo(# Bibliographical Information
        name="CASA",
        longName="Carnegie-Ames-Stanford approach", 
        version="1",
        entryAuthor="Verónika Ceballos-Núñez",
        entryAuthorOrcid="0000-0002-0046-1160",
        entryCreationDate="17/7/2015",
        doi="10.1029/93GB02725",
        further_references=BibInfo(doi="10.2307/1313568"),
        #  Also from PDF in Reflex experiment
        sym_dict=sym_dict
        
    ),
    #
    # the following variables constitute the compartmental system:
    #
    A,  # the overall compartmental matrix
    Input,  # the overall input
    t,  # time for the complete system
    x,  # state vector of the the complete system
    VegetationCarbonInputScalar(u),
    # vegetation carbon partitioning.
    VegetationCarbonInputPartitioningTuple(b),
    VegetationCarbonStateVariableTuple((C_f, C_w, C_r)),
    np1
})
