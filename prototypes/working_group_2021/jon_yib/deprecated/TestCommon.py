# Parent class to be used by different subclasses to use the same setUp
import sys
from pathlib import Path
import json 
import numpy as np
from testinfrastructure.InDirTest import InDirTest

sys.path.insert(0,'..')

from model_specific_helpers import (
    Observables,
    EstimatedParameters, 
    UnEstimatedParameters, 
    Parameters, 
    StateVariables,
    ModelParameters,
    get_example_site_vars,
    monthly_to_yearly,
)

config_file_name = 'config.json'
with Path(config_file_name).open(mode='r') as f:
    conf_dict=json.load(f) 
    dataPath = Path(conf_dict['dataPath'])
    
class TestCommon(InDirTest):
    def setUp(self):
        self.dataPath = dataPath
        npp, rh, ra, csoil, cveg = get_example_site_vars(Path(conf_dict['dataPath']))
        print(list(map(lambda a:a.shape,(npp,rh,ra,csoil,cveg))))
        #nyears = 320
        nyears = 2
        obs_tup=Observables(
            c_veg=cveg,
            c_soil=csoil,
            a_respiration=monthly_to_yearly(ra),
            h_respiration=monthly_to_yearly(rh)
        )

        self.obs = np.stack(obs_tup, axis=1)[0:nyears,:]
        
        self.epa0 = EstimatedParameters(
            beta_leaf=0.15,             # 0         
            beta_root=0.2,              # 1      
            k_leaf=1/365,               # 2      
            k_root=1/(365*5),           # 3         
            k_wood=1/(365*40),          # 4
            k_cwd=1/(365*5),            # 5      
            k_samet=0.5/(365*0.1),      # 6      
            k_sastr=0.5/(365*0.1),      # 7      
            k_samic=0.3/(365*0.137),    # 8      
            k_slmet=0.3/(365),          # 9      
            k_slstr=0.3/(365),          # 10      
            k_slmic=0.3/(365),          # 11      
            k_slow=0.3/(365*5),         # 12      
            k_arm=0.3/(222*365),        # 13      
            f_samet_leaf=0.3,           # 14      
            f_slmet_root=0.3,           # 15      
            f_samic_cwd=0.3,            # 16     
            C_leaf_0=cveg[0]/5,         # 17      
            C_root_0=cveg[0]/5,         # 18      
            C_cwd_0=cveg[0]/50,         # 19      
            C_samet_0=cveg[0]/300,      # 20      
            C_sastr_0=cveg[0]/300,      # 21      
            C_samic_0=cveg[0]/500,      # 22      
            C_slmet_0=csoil[0]/10,      # 23      
            C_slstr_0=csoil[0]/10,      # 24      
            C_slmic_0=csoil[0]/10,      # 25      
            C_slow_0=csoil[0]/10        # 26 
        )

        self.cpa = UnEstimatedParameters(
            C_soil_0=csoil[0],
            C_veg_0=cveg[0],
            rh_0 = rh[0],
            ra_0 = ra[0],
            npp=npp,
            clay=0.2028,
            silt=0.2808,
            nyears=nyears
        )
        self.pa = Parameters.from_EstimatedParametersAndUnEstimatedParameters(self.epa0,self.cpa)

        self.mpa = ModelParameters(**{k:v for k,v in self.pa._asdict().items() if k in ModelParameters._fields})
        
        self.epa0 = self.epa0 

