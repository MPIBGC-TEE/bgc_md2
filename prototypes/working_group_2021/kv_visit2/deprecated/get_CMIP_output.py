#!/usr/bin/env python3

# Description: script for downloading the CMIP data without using the user interface.  Requires user to set the
# source_id (model), variable_id (carbon pools), and data node (if the node set does not work).
# to get correct variable and model names check in the search filters at https://esgf-node.llnl.gov/search/cmip6/

# written by Alison C Bennett 30/11/2021 with code adapted from https://claut.gitlab.io/man_ccia/lab2.html and
# Adapted from: https://stackoverflow.com/a/37573701

## import packages
from pyesgf.search import SearchConnection # install using conda with: conda install -c conda-forge esgf-pyclient
import requests
import json
import os
from pathlib import Path
from tqdm import tqdm
from functools import reduce

## initialise connection (change node if desired)
conn = SearchConnection('https://esgf-node.llnl.gov/esg-search', distrib=True)

## set output directory
with Path('../config.json').open(mode='r') as f:
    conf_dict=json.load(f)
dp = Path(conf_dict["dataPath"])
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
