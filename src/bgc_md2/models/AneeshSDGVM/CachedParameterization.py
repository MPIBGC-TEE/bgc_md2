
import numpy as np
from ... import helper as h
from ...VirtualCachedParameterization import VirtualCachedParameterization
from collections import namedtuple

Drivers=namedtuple(
    "Drivers",
    ["npp"]
)    

def make_func_dict(dvs, **kwargs):
    return {
        "NPP": h.make_interpol_of_t_in_days(dvs.npp),
        "xi": lambda t: 1
    }

class CachedParameterization(VirtualCachedParameterization):
    Drivers = Drivers

    @property
    def func_dict(self):
        return make_func_dict(self.drivers)
