# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.13.6
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# make changes to imported files immidiately available 
# avoiding the need to reload (in most cases)
# %load_ext autoreload
# %autoreload 2
import general_helpers as gh
from pathlib import Path
variables=["cVeg", "cLitter", "cSoil","gpp","npp","ra","rh","tas","mrso","rh_annual","tsl","lai"]
models=np.array(("CLASSIC","CLM5","DLEM","IBIS","ISAM","ISBA_CTRIP",
"JSBACH","JULES-ES-1.0","LPJ-GUESS","LPJwsl","LPX-Bern","OCN",
"ORCHIDEE","ORCHIDEE-CNP","ORCHIDEEv3","ORCHIDEEv3_0.5deg",
"SDGVM","VISIT","YIBs"))
for i in range(len(models)):
    dataPath="C:/Users/konst/OneDrive - Cornell University/Data/Matrix MIP data/TRENDY/"+models[i]
    print(models[i])
    print(dataPath)
    gh.download_TRENDY_output(
        username="trendy-v9",
        password="gcb-2020",
        dataPath=Path(dataPath),#platform independent path desc. (Windows vs. linux)
        models=[models[i]],
        variables = variables, #Observables._fields + OrgDrivers._fields
        experiments = ["S2", "S3"]
    )



