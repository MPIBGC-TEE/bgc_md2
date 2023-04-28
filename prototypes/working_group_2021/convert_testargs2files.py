import general_helpers as gh
model_folders = [
    # first tier (best shape)
    "kv_visit2",  
    "jon_yib",
    "yz_jules",
    ##
    "Aneesh_SDGVM", #second tier (not quite ready)
    #"kv_ft_dlem",
    ##
    ##third tier
    ##"cj_isam", # has problems with ODE solution probably due to wrong parameters
    ## msh.numericX0 also yields a negative pool value for the last pool
    #"bian_ibis2",# 
    ##"cable-pop", 
]
def write_parameterization_from_test_args(
        test_args, 
        func_dict_param_dict,
        CachedParameterization: type, # a model specific class with a common interface
        data_path
    ):
    mvs = test_args.mvs
    svt = mvs.get_StateVariableTuple()
    V_init=test_args.V_init
    X_0_dict={k: V_init.__getattribute__(str(k)) for k in svt}
    cp = CachedParameterization(
        parameter_dict=test_args.par_dict,
        drivers=test_args.dvs,
        X_0_dict=X_0_dict,
        func_dict_param_dict=func_dict_param_dict
    )
    cp.write(data_path)

def write_data_assimilation_parameters_from_test_args_single(
        test_args, 
        data_path
    ):
    data_path.mkdir(exist_ok=True, parents=True)
    for s in ['cpa', 'epa_0', 'epa_min', 'epa_max', 'epa_opt']:
        h.dump_dict_to_json_path(
            (test_args.__getattribute__(s))._asdict(),
            data_path.joinpath(f"{s}.json"),
            indent=2
        )

def write_data_assimilation_parameters_from_test_args():
    for mf in model_folders:
        msh = gh.msh(mf)
        mvs=import_module(f"{msh.model_mod}.source").mvs
        CachedParameterization = import_module(f"{msh.model_mod}.CachedParameterization").CachedParameterization
        dir_path = Path(mf).joinpath("data_assimilation_parameters_from_test_args")
        try:
            cp=CachedParameterization.from_path(dir_path)
        except:    
            test_args = gh.test_args_2(mf)

            write_data_assimilation_parameters_from_test_args_single(
                test_args,
                CachedParameterization,
                dir_path
            )

