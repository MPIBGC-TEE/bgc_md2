
import netCDF4 as nc
import numpy as np
from ... import helper as h
from ...VirtualCachedParameterization import VirtualCachedParameterization
from collections import namedtuple

# Driver data streams on TRENDY server
Drivers = namedtuple(
    "Drivers",
    ["npp", "mrsos", "tsl", "landCoverFrac"]
)


def xi(tsl, Mw, Ms, Topt, Tcons, mrsos, landCoverFrac):
    # Notes:
    # 1.) fixme mm 4-11-2023: this function could be formulated as a
    # sympy piecewise expression. The parameters Mw, Ms, Topt,
    # Tcons,would then be part of the parameter_dict as all the others
    # and we would not need to treat them in a special way.
    # This would be a
    # bit easier to read since xi would be visible in the symbolic
    # representation.  Here we treat the xi parameters as an extra set

    # 2.)alternative FT
        # Q10 function (this is not what Clark et al 2011 Fig. 2 presented, the equation must be wrong)
    #FT = 2.0 ** ((tsl[mi] - 298.15) / 10)  # temperature rate modifier
        # RothC temperature function (Jenkinson 1990)
    FT = Tcons / (1 + np.exp(106/(tsl - 273.1 + Topt)))
    FV = 0.6 + 0.4 * (1 - landCoverFrac / 100)  # effect of vegetation cover
    # Mw is soil moisture at wilting point as a fraction of saturation
    # Ms is soil moisture content at saturation
    S0 = 0.5 * (1 + Mw)  # optimum soil moisture
    Smin = 1.7 * Mw  # lower threshold soil moisture for soil respiration
    if S0 < mrsos/Ms:
        FS = 1 - 0.8 * (mrsos/Ms - S0)  # effect of soil moisture
    if (Smin < mrsos/Ms) and (mrsos/Ms <= S0):
        FS = 0.2 + 0.8 * (mrsos/Ms - Smin) / (S0 - Smin)
    if mrsos/Ms <= Smin:
        FS = 0.2
    # print("FT,FV,FS", FT, FV, FS)
    rh_factor = FT * FV * FS
    return rh_factor # 1.0     # Set to 1 if no scaler implemented
    # return 1.0


def make_func_dict_2(dvs, Mw, Ms, Topt, Tcons):
        tsl_f, mrso_f, landCoverFrac_f, npp_f = map(
            h.make_interpol_of_t_in_days,
            (dvs.tsl, dvs.mrsos, dvs.landCoverFrac, dvs.npp)
        )
        return {
            "NPP": npp_f,
            "xi": lambda t: xi(
                tsl_f(t),
                Mw,
                Ms,
                Topt,
                Tcons,
                mrso_f(t),
                landCoverFrac_f(t)
             )
        }


class CachedParameterization(VirtualCachedParameterization):
    Drivers = Drivers
    func_dict_param_keys=("Mw","Ms","Topt","Tcons")
    @property
    def func_dict(self):
        fdps = self.func_dict_param_dict
        # add special keyword parameters
        d = make_func_dict_2(self.drivers, **fdps)
        return d
