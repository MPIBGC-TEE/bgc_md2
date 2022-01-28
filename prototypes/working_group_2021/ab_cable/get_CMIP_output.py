#!/usr/bin/env python3

# Description: script for downloading the CMIP data without using the user interface.  Requires user to set the
# source_id (model), variable_id (carbon pools), and data node (if the node set does not work).
# to get correct variable and model names check in the search filters at https://esgf-node.llnl.gov/search/cmip6/

# written by Alison C Bennett 30/11/2021 with code adapted from https://claut.gitlab.io/man_ccia/lab2.html and
# Adapted from: https://stackoverflow.com/a/37573701

## import packages
from pyesgf.search import SearchConnection # install using conda with: conda install -c conda-forge esgf-pyclient
import os
import pandas as pd
import requests
import json
from tqdm import tqdm

## initialise connection (change node if desired)
conn = SearchConnection('https://esgf-node.llnl.gov/esg-search', distrib=True)

## set output directory
os.chdir(dataPath)  # need to check this works with the config.json

## Set search query (model specific)
query= conn.new_context(
    latest= True, # don't change - gets latest version
    #facets = 'null',
    project='CMIP6', # don't change
    experiment_id='1pctCO2', # don't change
    source_id='ACCESS-ESM1-5', # Fixme change to your model name here
    variable_id='cCwd, cLeaf, cRoot', # Fixme set carbon pool names here
    data_node='esgf.nci.org.au')  # set the data node here - otherwise we get results for all datanodes and need to filter later.

n_files = query.hit_count

## get filename and url for all files

#initialise list
files = []

#loop through the files
for i in range(0, n_files):
    print('file ', i)
    hit = results[i].file_context().search()
    files2 = map(lambda f: {'filename': f.filename, 'url': f.download_url}, hit)
    files2 = list(files2)
    files2[0]
    files.extend(files2)

# extract to dataframe
files = pd.DataFrame.from_dict(files)

## download the data

#define the download function
def download(url, filename):
    print("Downloading ", filename)
    r = requests.get(url, stream=True)
    total_size, block_size = int(r.headers.get('content-length', 0)), 1024
    with open(filename, 'wb') as f:
        for data in tqdm(r.iter_content(block_size),
                         total=total_size // block_size,
                         unit='KiB', unit_scale=True):
            f.write(data)

    if total_size != 0 and os.path.getsize(filename) != total_size:
        print("Downloaded size does not match expected size!\n",
              "FYI, the status code was ", r.status_code)

# download the data
for index, row in files.iterrows():
    if os.path.isfile(row.filename):
        print("File exists. Skipping.")
    else:
        download(row.url, row.filename)