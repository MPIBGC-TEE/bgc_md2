import numpy as np
from ... import helper as h
from ...VirtualCachedParameterization import VirtualCachedParameterization
from collections import namedtuple

Drivers=namedtuple(
    "Drivers",
    ["npp","tas","gpp"]
)
        
def make_func_dict(dvs, **kwargs):

    def xi_leaf(tas):
        t_ref = 273.15 + 24
        t_half = 273.15 + 33
        t_exp = 1.8
        tf_frac = 0.2
        s_t = t_exp ** ((tas - t_ref)/10)
        s_f = (1 + np.exp(tf_frac * (tas-t_half)))
        return s_t / s_f 

    def xi_soil(tas):
        t_ref = 273.15 + 28
        t_half = 273.15 + 0
        t_exp = 1.9
        s_t = t_exp ** ((tas - t_ref)/10)
        s_f = 1 / (1 + np.exp(t_half - tas))
        return s_t * s_f 

    gpp_func, npp_func, tas_func = map(
        h.make_interpol_of_t_in_days,
        (dvs.gpp, dvs.npp, dvs.tas)
    )

    return {
        "temp": tas_func,
        "GPP": gpp_func,
        "NPP": npp_func,
        "xi_leaf": lambda t: xi_leaf(tas_func(t)),
        "xi_soil": lambda t: xi_soil(tas_func(t))
    }


class CachedParameterization(VirtualCachedParameterization):
    Drivers = Drivers

    @property
    def func_dict(self):
        return make_func_dict(self.drivers)
