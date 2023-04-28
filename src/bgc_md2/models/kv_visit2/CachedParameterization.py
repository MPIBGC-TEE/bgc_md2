from collections import namedtuple
from ... import helper as h
from ...VirtualCachedParameterization import VirtualCachedParameterization


OrgDrivers = namedtuple("OrgDrivers", ["gpp", "ra", "mrso", "tas"])
Drivers = namedtuple("Drivers", ("npp",) + OrgDrivers._fields[2:])

def make_func_dict(dvs,**kwargs):
    return {
        "TAS": h.make_interpol_of_t_in_days(dvs.tas),
        "mrso": h.make_interpol_of_t_in_days(dvs.mrso),
        "NPP": h.make_interpol_of_t_in_days(dvs.npp),
    }

class CachedParameterization(VirtualCachedParameterization):
    Drivers = Drivers

    @property
    def func_dict(self):
        return make_func_dict(self.drivers)
