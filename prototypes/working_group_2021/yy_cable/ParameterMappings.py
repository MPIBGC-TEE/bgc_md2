from collections import namedtuple
UnEstimatedParameters = namedtuple(
    "UnEstimatedParameters",
    [
        'C_leaf_0',
        'C_root_0',
        'C_wood_0',
        'clitter_0',
        'csoil_0',
        'rh_0',
        'clay',
        'silt',
        'lig_wood',
        'f_wood2CWD',
        'f_metlit2mic',
        'npp',
        'number_of_months'
    ]
)

# This set is used (as argument) by the functions that are called
# inside the mcmc
EstimatedParameters = namedtuple(
    "EstiamatedParameters",
    [
        "beta_leaf",    #  0 (indices uses in original code) 
        "beta_root",    #  1
        "lig_leaf",     #  2
        "f_leaf2metlit",#  3
        "f_root2metlit",#  4
        "k_leaf",       #  5
        "k_root",       #  6
        "k_wood",       #  7
        "k_metlit",	#  8
        "k_mic",	#  9
        "k_slowsom",	# 10
        "k_passsom",	# 11
        "C_metlit_0",	# 12
        "C_CWD_0",	# 13
        "C_mic_0",	# 14
        "C_passom_0"    # 15
    ]
)

# This is the set off all 
_Parameters=namedtuple(
        "Parameters",
        EstimatedParameters._fields + UnEstimatedParameters._fields
)
class Parameters(_Parameters):
    @classmethod
    def from_EstimatedParametersAndUnEstimatedParameters(
            cls,
            epa :EstimatedParameters,
            cpa :UnEstimatedParameters
        ):
        return cls(*(epa + cpa))


# This set defines the order of the c pools
# The order is crucial for the compatibility
# with the matrices (B and b) If you change ist
# the matrix changes
StateVariables = namedtuple(
    'StateVariables',
    [
        'C_leaf',
        'C_root',
        'C_wood',
        'C_metlit',
        'C_stlit',
        'C_CWD',
        'C_mic',
        'C_slowsom',
        'C_passsom'
    ]
)
Observables = namedtuple(
    'Observables',
    [
        'C_leaf',
        'C_root',
        'C_wood',
        'c_litter',
        'c_soil',
        'respiration'
    ]
)

# We define another set of parameters which describes
# the parameters of the matrices A,K and the vector b
# and drivers like npp (in form of arrays)
# but does not include start values and hyperparameters like the 'number_of_months'
# This distinction is helpful for the forward simulation where the
# distinction between estimated and constant is irrelevant.
ModelParameters = namedtuple(
    "ModelParameters",
    [
        name for name in Parameters._fields 
        if name not in [
            'C_leaf_0',
            'C_root_0',
            'C_wood_0',
            'clitter_0',
            'csoil_0',
            "C_metlit_0",
            "C_CWD_0",
            "C_mic_0",
            "C_passom_0",
            'number_of_months'
        ] 
    ]
)
