import numpy as np
from copy import copy
from sympy import Symbol, ImmutableMatrix, Function, sympify
import CompartmentalSystems.helpers_reservoir as hr
from CompartmentalSystems.smooth_reservoir_model import SmoothReservoirModel
from CompartmentalSystems.smooth_model_run import SmoothModelRun
import CompartmentalSystems.start_distributions as sd
from . import mvars 
def sub_mr_smav(
        smr: SmoothModelRun,
        #start_mean_age_vec: mvars.NumericStartMeanAgeTuple,
        svt
    ):
    srm = copy(smr.model)
    sv = srm.state_vector
    t = srm.time_symbol
    svl = list(sv)
    svs = set(sv)
    #from IPython import embed; embed()
    X_fix=smr.start_values
    times = smr.times
    ifbs = hr.in_fluxes_by_symbol(sv, srm.input_fluxes)
    intfbs= hr.internal_fluxes_by_symbol_2(sv,srm.internal_fluxes)
    ofbs= hr.out_fluxes_by_symbol_2(sv,srm.output_fluxes)
    combined = (
        set(sv),
        ifbs,
        ofbs,
        intfbs,
    )

    (
        sub_sv_set,
        sub_in_fluxes,
        sub_out_fluxes,
        sub_internal_fluxes,
    ) = hr.extract(combined, set(svt))

    # we have to provide solutions (functions of time)  for the
    # state variables of the outer system ( In our case soil
    # influxes depend on the poolvalues of the veg system. In
    # general every subsystem could depend on pools of every other
    # subsystem
    outer_pools = svs.difference(sub_sv_set)
    # first replace the symbols by functions
    subs_dict = {sym: Function(str(sym))(t) for sym in outer_pools}
    sub_in_fluxes_f, sub_internal_fluxes_f, sub_out_fluxes_f = map(
        lambda d: {k: sympify(flux_ex).subs(subs_dict) for k, flux_ex in d.items()},
        (sub_in_fluxes, sub_internal_fluxes, sub_out_fluxes),
    )

    # now prepare the extension of the func_dict by the solutions
    # for the variables

    #start_mean_age_vec2 = sd.start_age_moments_from_steady_state(
    #    smr.model,
    #    t0=np.min(times),
    #    parameter_dict=smr.parameter_dict,
    #    func_set=smr.func_set,
    #    max_order=1,
    #).reshape(-1)
    order = 0
    s_arr, s_func = smr._solve_age_moment_system(
        order,
        None#start_mean_age_vec.reshape(1, smr.nr_pools),
    )

    #from IPython import embed; embed()
    def s_func_maker(sym):
        return lambda t: s_func(t)[svl.index(sym)]

    outer_pool_funcs = {
        # note that we create the function in another
        # function since the normal dictionary comprehension
        # fails due to pythons lazy evaluation not
        # protecting the function specific variable
        # (in this case sym)
        #Function(str(sym)): s_func_maker(sym)
        str(sym): s_func_maker(sym)
        for sym in outer_pools
    }
    sub_func_dict = {**smr.func_set, **outer_pool_funcs}
    sub_X_fix = X_fix[[svl.index(w) for w in svt]]
    #sub_mvs = CMTVS(
    #    {
    #        t,
    #        StateVariableTuple(svt),
    #        InFluxesBySymbol(sub_in_fluxes_f),
    #        InternalFluxesBySymbol(sub_internal_fluxes_f),
    #        OutFluxesBySymbol(sub_out_fluxes_f),
    #        NumericParameterization(
    #            par_dict=par_dict,
    #            func_dict=sub_func_dict,
    #        ),
    #        NumericStartValueArray(sub_X_fix),
    #        NumericSimulationTimes(times),
    #    },
    #    mvs.computers,
    #)
    #sub_smr = sub_mvs.get_SmoothModelRun()
    sub_srm = SmoothReservoirModel.from_state_variable_indexed_fluxes(
        state_vector=svt,
        time_symbol=t,
        input_fluxes=sub_in_fluxes_f,
        output_fluxes=sub_out_fluxes_f,
        internal_fluxes=sub_internal_fluxes_f,
    )
    sub_smr = SmoothModelRun(
            sub_srm,
            smr.parameter_dict,
            sub_X_fix,
            times,
            sub_func_dict
    )
    # for the start mean ages we can not just take the respective parts of the system mean ages
    # since they refer to the times stuff has spent in the SYSTEM and not in the SUB system.
    # Although the pools in the subsystem are pools in the
    #sub_start_mean_age_vec = sd.start_age_moments_from_steady_state(
    #    sub_smr.model,
    #    t0=np.min(times),
    #    parameter_dict=smr.parameter_dict,
    #    func_set=sub_func_dict,
    #    max_order=1,
    #)
    return sub_smr 
