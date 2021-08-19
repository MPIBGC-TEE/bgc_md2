#!/usr/bin/env python
from mm_helpers import get_example_site_vars, mcmc
import numpy as np
import matplotlib 
from matplotlib import pyplot as plt

dataPath = Path('/home/data/yuanyuan/')
data = get_example_site_vars(dataPath)
nyears =140 
df, df_j = mcmc(dataPath,start_pa=pa1,nyears=140)
df.to_csv(str(dataPath.joinpath('cable_demo_da_aa.csv')),sep=',')
df_j.to_csv(str(dataPath.joinpath('cable_demo_da_j_aa.csv')),sep=',')
