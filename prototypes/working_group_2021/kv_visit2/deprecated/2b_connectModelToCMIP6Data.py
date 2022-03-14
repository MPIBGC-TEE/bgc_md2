# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light
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

# ## Connecting symbolic description and data
#
# ### Intermediate summary:
# We have achieved the symbolic formulation of the model. We can use it to check the structure and compute derived diagnostic variables including the matrices used for the traceability analysis, but up to now only in symbolic form. 
#
# The next goal is to connect the symbolic formulation to the data.
# Since the data comes in several shapes this involves several steps. 
# We also want to be able to make this step portable across different models and computers.
# The endproduct is a collection of models that everybody can run who installs the package and executes the code we provide
#
# We will have to: 
# 1. Find as many model parameters as possible in the model description (in the literature or in communication with the modeling group) so that we do not have to estimate them from the model output. 
# 1. provide code to download the output for your model.
# 1. implement functions for the drivers (using the data)
# 1. run the model forward with a possible set of parameters.
# 1. infer unknown parameters by data assimilation.
#
# ### downloading the data
# #### create a small site specific config file 
# This file specifies:
# - a username and password to download the data 
# - the location where you want to download the data to 
#   which will differ depending on the machine you are using (your laptop or a supercomputer) and     also accross users. You will have to have one everywhere you want to work with the model.
# Here comes a template from my laptop (content of file `../config.json`):
# `{"username": "trendy-v9", "password": "gcb-2020", "dataPath": "/home/data/VISIT_data_CMIP6"}`
# - Please **do not** put your `config.json` directly under version control,   since it is different for every user and will only create merge conflicts if you do. Rather copy it to `config.json.{yourname}` and add this file to git. So  everybody can copy and adapt it.
#
# Note that 
# - the file resides one level above your current folder since it is not modelspecific
#   (This is a change from the first round of model gathering)

# load the model as we have specified it before in source.py
from source import mvs
import bgc_md2.display_helpers as dh
import bgc_md2.helper as h

h.compartmental_graph(mvs)

dh.mass_balance_equation(mvs)

# +
import sys
sys.path.insert(0,'..') # necessary to import general_helpers
#from general_helpers import download_TRENDY_output
import json 
from pathlib import Path

with Path('../config.json').open(mode='r') as f:
    conf_dict=json.load(f) 
    
dp=Path(conf_dict['dataPath'])


# +
#this function will later be stored in a file model `specific_helpers.py`
# you have to adapt it to download the data for YOUR model

def download_my_CMIP6_data(dp):
    # Description: downloading the CMIP data without using the user interface.  Requires user to set the
    # source_id (model), variable_id (carbon pools), and data node (if the node set does not work).
    # to get correct variable and model names check in the search filters at https://esgf-node.llnl.gov/search/cmip6/
    
    # written by Alison C Bennett 30/11/2021 with code adapted from https://claut.gitlab.io/man_ccia/lab2.html and
    # Adapted from: https://stackoverflow.com/a/37573701
    
    from pyesgf.search import SearchConnection # install using conda with: conda install -c conda-forge esgf-pyclient
    import requests
    import json
    import os
    from pathlib import Path
    from tqdm import tqdm
    from functools import reduce
    
    ## initialise connection (change node if desired)
    conn = SearchConnection('https://esgf-node.llnl.gov/esg-search', distrib=True)
    
    if not dp.exists():
        dp.mkdir(exist_ok=True)
    
    ## Set search query (model specific)
    query = conn.new_context(
        latest= True, # don't change - gets latest version
        #facets = 'null',
        project='CMIP6', # don't change
        experiment_id='1pctCO2', # don't change
        source_id='ACCESS-ESM1-5', # Fixme change to your model name here
        variable_id="cLeaf, cLitterAbove, cLitterBelow, cRoot, cSoilFast, cSoilMedium, cSoilSlow, cVeg, fLitterSoil, fVegLitter, mrsos, npp, rh, tsl" ,
        data_node='esgf.nci.org.au')  # set the data node here - otherwise we get results for all datanodes and need to filter later.
    
    n_files = query.hit_count
    results=query.search()
    ## get filename and url for all files
    
    def dictlist(result):
        hit_set = result.file_context().search()
        return list(map(lambda f: {'filename': f.filename, 'url': f.download_url}, hit_set))
    
    file_dicts = reduce(
        lambda acc,r : acc + dictlist(r),
        results,
        []
    )
    #from IPython import embed; embed()
    ###define the download function
    def download(url, filename):
        print("Downloading ", filename)
        print("url", url)
        r = requests.get(url, stream=True)
        total_size, block_size = int(r.headers.get('content-length', 0)), 1024
        p=dp.joinpath(filename)
        with p.open('wb') as f:
            for data in tqdm(r.iter_content(block_size),
                             total=total_size // block_size,
                             unit='KiB', unit_scale=True):
                f.write(data)
    
        if total_size != 0 and os.path.getsize(p) != total_size:
            print("Downloaded size does not match expected size!\n",
                  "FYI, the status code was ", r.status_code)
    
    ### download the data
    for d in file_dicts:
        if dp.joinpath(d['filename']).exists():
            print("File exists. Skipping.")
        else:
            download(d["url"], d["filename"])
    

download_my_CMIP6_data(dp)

# +

import netCDF4 as nc
import numpy as np
from pathlib import Path
import json 
def get_variables_from_files(dataPath):
    # Read NetCDF data  ******************************************************************************************************************************
    names = [
        ('cLeaf', 'cLeaf_Lmon_MIROC-ES2L_1pctCO2-bgc_r1i1p1f2_gn_185001-199912.nc'),
        ('cLitterAbove', 'cLitterAbove_Lmon_MIROC-ES2L_1pctCO2-bgc_r1i1p1f2_gn_185001-199912.nc'),
        ('cLitterBelow', 'cLitterBelow_Lmon_MIROC-ES2L_1pctCO2-bgc_r1i1p1f2_gn_185001-199912.nc'),
        ('cRoot', 'cRoot_Lmon_MIROC-ES2L_1pctCO2-bgc_r1i1p1f2_gn_185001-199912.nc'),
        ('cSoilFast', 'cSoilFast_Lmon_MIROC-ES2L_1pctCO2-bgc_r1i1p1f2_gn_185001-199912.nc'),
        ('cSoilMedium', 'cSoilMedium_Lmon_MIROC-ES2L_1pctCO2-bgc_r1i1p1f2_gn_185001-199912.nc'),
        ('cSoilSlow', 'cSoilSlow_Lmon_MIROC-ES2L_1pctCO2-bgc_r1i1p1f2_gn_185001-199912.nc'),
        ('cVeg', 'cVeg_Lmon_MIROC-ES2L_1pctCO2-bgc_r1i1p1f2_gn_185001-199912.nc'),
        ('fLitterSoil', 'fLitterSoil_Lmon_MIROC-ES2L_1pctCO2-bgc_r1i1p1f2_gn_185001-199912.nc'),
        ('fVegLitter', 'fVegLitter_Lmon_MIROC-ES2L_1pctCO2-bgc_r1i1p1f2_gn_185001-199912.nc'),
        ('mrsos', 'mrsos_Lmon_MIROC-ES2L_1pctCO2-bgc_r1i1p1f2_gn_185001-199912.nc'),
        ('npp', 'npp_Lmon_MIROC-ES2L_1pctCO2-bgc_r1i1p1f2_gn_185001-199912.nc'),
        ('rh', 'rh_Lmon_MIROC-ES2L_1pctCO2-bgc_r1i1p1f2_gn_185001-199912.nc'),
        ('tsl', 'tsl_Lmon_MIROC-ES2L_1pctCO2-bgc_r1i1p1f2_gn_185001-199912.nc'),
    ]
    def f(tup):
        vn, fn = tup
        path = dataPath.joinpath(fn)
        ds = nc.Dataset(str(path))
        return ds.variables[vn][:, :, :]

    return map(f, names)

#     # Read NetCDF data  ******************************************************************************************************************************

def 
(dataPath):
    (
        C_leaf,
        C_litter_above,
        C_litter_below,
        C_root,
        C_fast_som,
        C_slow_som,
        C_pass_som,
        C_veg,
        f_litter2som,
        f_veg2litter,
        mrso,
        npp,
        rh,
        tsl
    )= get_variables_from_files(dataPath)
    # pick up 1 site   wombat state forest
    s = slice(None, None, None)  # this is the same as :
    t = s, 50, 33  # [t] = [:,49,325]
    npp = npp[t] * 86400   # kg/m2/s kg/m2/day;
    rh = rh[t]*86400  # per s to per day
    f_veg2litter = f_veg2litter[t] * 86400
    f_litter2som = f_litter2som[t] * 86400
    tsl_mean = np.mean(tsl, axis=1)  # average soil temperature at different depth
    (
        C_leaf,
        C_litter_above,
        C_litter_below,
        C_root,
        C_fast_som,
        C_slow_som,
        C_pass_som,
        C_veg,
        mrso,
        tsl
    ) = map(
        lambda var: var[t],
        (
            C_leaf,
            C_litter_above,
            C_litter_below,
            C_root,
            C_fast_som,
            C_slow_som,
            C_pass_som,
            C_veg,
            mrso,
            tsl_mean
        )
    )
    C_wood = C_veg - C_leaf - C_root
    return (npp, C_leaf, C_wood, C_root, C_litter_above, C_litter_below, C_fast_som, C_slow_som, C_pass_som,
            rh, f_veg2litter, f_litter2som, mrso, tsl)


with Path('config.json').open(mode='r') as f:
    conf_dict=json.load(f) 

dataPath = Path(conf_dict['dataPath'])
(
    npp,
    C_leaf,
    C_wood,
    C_root,
    C_litter_above,
    C_litter_below,
    C_fast_som,
    C_slow_som,
    C_pass_som,
    rh,
    f_veg2litter,
    f_litter2som,
    mrso,
    tsl
)=get_example_site_vars(dataPath)

import sys 
sys.path.insert(0,'..')
from general_helpers import day_2_month_index
def NPP_fun(day ):
    return npp[day_2_month_index(day)] 

func_dict={NPP: NPP_fun}
# -

# ### Forward run
# The next goal is to run the model forward with a given set of parameters.
# So we need:
# 1. values for all the parameters 
# 1. implementations for the symbolic functions 
# 1. start values for all the pool contents

# In this example we have the initial values for the elements of the K and A matrices ($k_i$ and $f_{i,j}$ ) but we want the values for the individual flux rates. 
# So we first create the symbolic $k_{i}$ and $f_{i,j}$ which gives us an alternative description 
# of the product $M = K A=K_{sym} A_{sym}$ from this matrix we can compute the outfluxes and internal fluxes. If we assume that the same fluxes could be expressed as $flux=fluxrate * donorpoolcontent$ we can compute the fluxrates 
#
# **This step is hopefully not necessarry for your model since you should find parameters directly**
# Since we unfortunately started with a different description, we have to convert 
#

# +
from sympy import ImmutableMatrix
sv=mvs.get_StateVariableTuple()
n=len(sv)
# create new symbols for the f_{i,j}
for i in range(n):
    for j in range(n):
        if A[i,j]!=0 and i!=j:
            name="f_" + str(sv[j]) + "_2_" + str(sv[i])
            code="{0}=Symbol('{0}')".format(name)
            print(code)
            exec(code)
            
A_sym=ImmutableMatrix(
    n,n,
    lambda i,j:  -1 if i==j else (
        0 if A[i,j]==0 else Symbol("f_" + str(sv[j]) + "_2_" + str(sv[i]))
    )
)
A_sym

# +
# create new symbols for the k_{i}
for i in range(n):
    if K[i,i]!=0:
        name="k_{0}".format(sv[i])
        code="{0}=Symbol('{0}')".format(name)
        print(code)
        exec(code)
        
K_sym=ImmutableMatrix(
    n,n,
    lambda i,j: Symbol("k_" + str(sv[i])) if i==j else 0
)
K_sym
# -

M_sym=A_sym*K_sym
M_sym

import  CompartmentalSystems.helpers_reservoir as hr
hr.out_fluxes_by_symbol(sv,M_sym)

# we create a dictionary for the outfluxrates (flux divided by dono pool content)
outflux_rates = {"k_"+str(key)+"_rh":value/key for key,value in hr.out_fluxes_by_symbol(sv,M_sym).items()}
internal_flux_rates = {"k_"+str(key[0])+"_2_"+str(key[1]):value/key[0] for key,value in hr.internal_fluxes_by_symbol(sv,M_sym).items()}
from copy import  deepcopy
all_rates=deepcopy(outflux_rates)
all_rates.update(internal_flux_rates)
all_rates

# +
# and one for the internal fluxrates
# -

old_par_dict = {
    beta_leaf: 0.6,
    beta_wood: 0.25,
    T_0: 2,
    E: 4,
    KM: 10,
    f_C_leaf_2_C_leaf_litter: 1,
    f_C_wood_2_C_wood_litter: 1,
    f_C_root_2_C_root_litter: 1,
    f_C_leaf_litter_2_C_soil_fast: 0.41,
    f_C_leaf_litter_2_C_soil_slow: 0.07,
    f_C_leaf_litter_2_C_soil_passive: 0.02,
    f_C_wood_litter_2_C_soil_fast: 0.30,
    f_C_wood_litter_2_C_soil_slow: 0.12,
    f_C_wood_litter_2_C_soil_passive: 0.08,
    f_C_root_litter_2_C_soil_fast: 0.30,
    f_C_root_litter_2_C_soil_slow: 0.14,
    f_C_root_litter_2_C_soil_passive: 0.07,
    k_C_leaf: 1 / (60 * 2),
    k_C_wood: 1 / (365 * 30),
    k_C_root: 1 / (365 * 22),
    k_C_leaf_litter: 1 / (365 * 3.3),
    k_C_wood_litter: 1 / (365 * 11),
    k_C_root_litter: 1 / (365 * 11),
    k_C_soil_fast: 1 / (365 * 18),
    k_C_soil_slow: 1 / (365 * 100),
    k_C_soil_passive: 1 / (365 * 350),
}


# Now we can translate the old paramterisation to the new one.

par_dict={
    beta_leaf: 0.6,
    beta_wood: 0.25,
    T_0: 2, 
    E: 4, 
    KM: 10 
}
par_dict.update(
    {k:v.subs(old_par_dict) for k,v in all_rates.items()}
)
par_dict

# To be able to run the model forward we not only have to replace parameter symbols by values but symbolic functions by normal python functions.
# In our case the functions for $NPP$ and $\xi$ have to be provided. NPP_fun will interpolate the NPP for the day in question from the data. Which we have to load. 
# We will later store these functions in  `model_specific_helpers.py` which resides in the same folder as this notebook. You will have to adapt them to your data set. 



To make a 


