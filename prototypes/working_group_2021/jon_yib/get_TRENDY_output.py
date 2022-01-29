#!/usr/bin/env python3

import paramiko
from pathlib import Path
import json
import tarfile

# open a transport
host = "trendy.ex.ac.uk"
port = 22
transport = paramiko.Transport(host)

with Path('config.json').open(mode='r') as f:
    conf_dict=json.load(f) 
# authentication
transport.connect(None,username=conf_dict["username"],password=conf_dict["password"])
sftp = paramiko.SFTPClient.from_transport(transport)

# download files
remote_path = "output"

model      = "YIBs"
experiment = "S2"
variables   = ["Monthly_npp","Monthly_rh", "Monthly_ra", "Annual_cSoil", "Annual_cVeg"]
ext = ".nc"
full_ext = ".tar.gz"

# create the directory for the data
Path(conf_dict['dataPath']).mkdir(exist_ok=True)
for variable in  variables:
    trunk = model + "_" + experiment + "_" + variable 
    filename  = trunk + ext + full_ext
    localpath = conf_dict['dataPath'] + "/" + filename
    try:
        sftp.get(
            remotepath=remote_path + "/" + model+ "/" + experiment + "/" + filename,
            localpath=localpath
        )
        f=tarfile.open(localpath,'r:gz')
        f.extractall(path=conf_dict['dataPath'])
        f.close()


    except FileNotFoundError:
        print("file not found:" + filename)

                    
print("finished!")
