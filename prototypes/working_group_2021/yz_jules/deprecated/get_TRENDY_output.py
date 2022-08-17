#!/usr/bin/env python3

import paramiko
from pathlib import Path
import json


# open a transport
host = "trendy.ex.ac.uk"
port = 22
transport = paramiko.Transport(host)


with Path('config.json').open(mode='r') as f:
    conf_dict = json.load(f)
# authentication
transport.connect(None, username=conf_dict["username"], password=conf_dict["password"])


sftp = paramiko.SFTPClient.from_transport(transport)

# download files
remote_path = "output"

# FIXME please update to refer to data location in config.json instead of hardcoded path here.
local_path = conf_dict["dataPath"]


# Other models, "CLASSIC","CLM5","DLEM","IBIS","ISAM","ISBA_CTRIP","JSBACH","JULES-ES","LPJ-GUESS","LPJwsl","LPX-Bern",
#                 "OCN","ORCHIDEEv3","SDGVM","VISIT","YIBs"

models      = ["JULES-ES"]
experiments = ["S2"]
variables   = ["cCwd", "cLeaf", "cLitter", "cRoot", "cSoil", "cVeg", "cWood", "npp", "rh"]

for model in models:
    print("downloading model ", model)
    for experiment in experiments:
        for variable in variables:

            modelname = model
            modelname_file = model
            ext = "nc"
            extra = ""
            
            if model == "CLM5":
                modelname_file = "CLM5.0"
            elif model == "ISBA_CTRIP":
                modelname_file = "ISBA-CTRIP"
            elif model == "JULES-ES":
                modelname = "JULES-ES-1.0"
                modelname_file = "JULES-ES-1p0"
            elif model == "SDGVM" or model == "VISIT":
                ext = "nc.gz"
            elif model == "YIBs":
                ext = "nc.tar.gz"
                extra = "Monthly_"
            elif model == "LPJwsl":
                modelname_file = "LPJ"
                ext = "nc.gz"
                
            filename  = modelname_file + "_" + experiment + "_" + extra + variable + "." + ext
               
            try:
                #print(sftp.stat(remote_path + "/" + modelname + "/" + experiment + "/" + filename)) # get file information + test if existing (IOError handling)
                sftp.get(remotepath=remote_path + "/" + modelname + "/" + experiment + "/" + filename,
                         localpath=local_path + "/" + filename)
            #except IOError:
            #    print("file not found:" + filename)
            except FileNotFoundError:
                print("file not found:" + filename)

                    
print("finished!")
