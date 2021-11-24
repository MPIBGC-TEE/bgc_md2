# Parent class to be used by different subclasses to use the same setUp
import sys
from pathlib import Path
import json 
from testinfrastructure.InDirTest import InDirTest

sys.path.insert(0,'..')

from model_specific_helpers import (
    EstimatedParameters, 
    UnEstimatedParameters, 
    Parameters, 
    StateVariables,
    ModelParameters,
    get_example_site_vars,
)

config_file_name = 'config.json'
with Path(config_file_name).open(mode='r') as f:
    conf_dict=json.load(f) 
    dataPath = Path(conf_dict['dataPath'])
    
class TestCommon(InDirTest):
    def setUp(self):
        self.dataPath = dataPath
        ####################################################
        # The estimated parameters include start values
        # for some of the pools
        epa0 = EstimatedParameters(
            beta_leaf=0.15,
	    beta_root=0.2,
	    lig_leaf=0.15,
	    f_leaf2metlit=0.28,
	    f_root2metlit=0.6,
	    k_leaf=1/365,
	    k_root=1/(365*5),
	    k_wood=1/(365*40),
	    k_metlit=0.5/(365*0.1),
	    k_mic=0.3/(365*0.137),
	    k_slowsom=0.3/(365*5),
	    k_passsom=0.3/(222.22*365),
	    C_metlit_0=0.05,
	    C_CWD_0=0.1,
	    C_mic_0=1,
	    C_passom_0=5,
        )


        npp, rh, clitter, csoil, cveg, cleaf, croot, cwood = get_example_site_vars(self.dataPath)

        self.x_init = StateVariables(
            C_leaf=cleaf[0],
            C_root=croot[0],
            C_wood=cwood[0],
            C_metlit=epa0[12],
            C_stlit=epa0[13],
            C_CWD=clitter[0]-epa0[12]-epa0[13],
            C_mic=epa0[14],
            C_slowsom=csoil[0]- epa0[14] - epa0[15],
            C_passsom=epa0[15]
        )
       
        cpa= UnEstimatedParameters(
            C_leaf_0 = cleaf[0],
            C_root_0 = croot[0],
            C_wood_0 = cwood[0],
            clitter_0 = clitter[0],
            csoil_0 = csoil[0],
            rh_0 = rh[0],
            clay=0.2028,
            silt=0.2808,
            lig_wood=0.4,
            f_wood2CWD=1,
            f_metlit2mic=0.45,
            npp=npp, #array!!!
            number_of_months=10
        )
        pa = Parameters.from_EstimatedParametersAndUnEstimatedParameters(epa0,cpa)

        mpa = ModelParameters(**{k:v for k,v in pa._asdict().items() if k in ModelParameters._fields})
        
        self.npp= npp
        self.epa0 = epa0 
        self.cpa = cpa
        self.pa = pa 
        self.mpa = mpa 
