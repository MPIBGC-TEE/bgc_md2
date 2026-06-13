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
            beta_leaf=0.6,  # 0 (parameters used in original code)
            beta_wood=0.25,  # 1
            k_leaf_litter_2_soil_fast=0.41 / (365 * 3.3),  # 2
            k_leaf_litter_2_soil_slow=0.07 / (365 * 3.3),  # 3
            k_leaf_litter_2_soil_passive=0.02 / (365 * 3.3),  # 4
            k_wood_litter_2_soil_fast=0.30 / (365 * 11),  # 5
            k_wood_litter_2_soil_slow=0.12 / (365 * 11),  # 6
            k_wood_litter_2_soil_passive=0.08 / (365 * 11),  # 7
            k_root_litter_2_soil_fast=0.30 / (365 * 11),  # 8
            k_root_litter_2_soil_slow=0.14 / (365 * 11),  # 9
            k_root_litter_2_soil_passive=0.07 / (365 * 11),  # 10
            k_leaf_2_leaf_litter =  1 / (60 * 2)   ,  # 11
            k_wood_2_wood_litter =  1 / (365 * 30) ,  # 12
            k_root_2_root_litter =  1 / (365 * 22) ,  # 13
            k_leaf_litter_rh= 1 - (.41+.07+.02)  / (365 * 3.3),  # 14
            k_wood_litter_rh= 1- (.3+.12+.08) / (365 * 11),  # 15
            k_root_litter_rh= 1- (.3+.14+.07) / (365 * 11),  # 16
            k_soil_fast_rh = 1 / (365 * 18),  # 17
            k_soil_slow_rh = 1 / (365 * 100),  # 18
            k_soil_passive_rh =1 / (365 * 350),  # 19
            C_leaf_lit_0=0.3,  # 20
            T_0=2,  # 21
            E=4,  # 22
            KM=10  # 23
        )


        nyears = 150  # reduced time span for testing purposes
        tot_len = 12 * nyears
        
        npp, C_leaf, C_wood, C_root, C_litter_above, C_litter_below, C_fast_som, C_slow_som, C_pass_som, \
        rh, F_veg2litter, F_litter2som, mrso, tsl = get_example_site_vars(dataPath)
        cpa = UnEstimatedParameters(
            C_leaf_0=C_leaf[0],
            C_wood_0=C_wood[0],
            C_root_0=C_root[0],
            C_litter_above_0=C_litter_above[0],
            C_litter_below_0=C_litter_below[0],
            C_soil_fast_0=C_fast_som[0],
            C_soil_slow_0=C_slow_som[0],
            C_soil_passive_0=C_pass_som[0],
            rh_0=rh[0],
            F_veg_lit_0=F_veg2litter[0],
            F_lit_soil_0=F_litter2som[0],
            npp=npp,
            number_of_months=tot_len,
            mrso=mrso,
            tsl=tsl
        )

        self.x_init = StateVariables(
            C_leaf=cpa.C_leaf_0,
            C_wood=cpa.C_wood_0,
            C_root=cpa.C_root_0,
            C_leaf_litter=epa0.C_leaf_lit_0,
            C_wood_litter=cpa.C_litter_above_0 - epa0.C_leaf_lit_0,
            C_root_litter=cpa.C_litter_below_0,
            C_soil_fast=cpa.C_soil_fast_0,
            C_soil_slow=cpa.C_soil_slow_0,
            C_soil_passive=cpa.C_soil_passive_0
        )
       
        pa = Parameters.from_EstimatedParametersAndUnEstimatedParameters(epa0,cpa)

        mpa = ModelParameters(**{k:v for k,v in pa._asdict().items() if k in ModelParameters._fields})
        
        self.npp= npp
        self.epa0 = epa0 
        self.cpa = cpa
        self.pa = pa 
        self.mpa = mpa 
