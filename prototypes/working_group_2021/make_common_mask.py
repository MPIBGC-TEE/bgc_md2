#!/usr/bin/env python
# %load_ext autoreload
# %autoreload 2
import numpy as np
import general_helpers as gh
import matplotlib.pyplot as plt
from pathlib import Path


# model_folders = [
#     'yz_jules',
#     'cj_isam',
#     'Aneesh_SDGVM',
#     'kv_ft_dlem',
#     'kv_visit2',
#     'jon_yib',
#     'bian_ibis2'
# ]
model_folders = [
                #"bian_ibis2", # to get a mask with deserts exclude IBIS models
                "ab_classic",
                "jsbach",
                "clm5",
                "ORCHIDEE-CNP","ORCHIDEEv3",
                "kv_ft_dlem","cj_isam","isba-ctrip",
                "yz_jules","lpj-guess","lpjwsl","lpx-bern",
                "ORCHIDEE-V2","Aneesh_SDGVM","kv_visit2","jon_yib",    
                "ORCHIDEE",
                #"ORCHIDEEv3_0.5deg", 
                #"CABLE_POP"          
                ]


def mask(mf):
    print("###########################################")
    print("#######   {}   ######".format(mf))
    msh = gh.msh(mf)
    return msh.spatial_mask(Path(gh.confDict(mf)['dataPath']))

# +
masks = list(
    map(
        mask,
        model_folders
    )
)
ma = np.zeros(
    shape=(360, 720),
    # dtype=np.bool_
)
ut = gh.CoordMask(
    index_mask=ma,
    tr=gh.globalMaskTransformers(ma)
)
masks.append(ut)
###############################################################################
# plotting
titles = model_folders+["empty_mask"]
n = len(masks)
f = plt.figure(figsize=(20,n*10))
for i in range(n):
    m = masks[i]
    ax = f.add_subplot(n+1, 1, i+1)
    ax.set_title(titles[i])
    m.plot_dots(ax)
    
# At the end of the list add the  ultimate target mask (empty)
common_mask = gh.combine_masks_2(masks)
# common_mask=gh.project_2(source=masks[0],target=masks[1])
# common_mask=gh.project_2(source=masks[0],target=ut)

#common_mask.write_netCDF4(Path("common_mask.nc"))
#common_mask.write_netCDF4(Path("common_mask_all_models.nc"))
common_mask.write_netCDF4(Path("common_mask_all_models_w_deserts.nc"))

ax = f.add_subplot(n+1, 1, n+1)
common_mask.plot_dots(ax)

# +
#f.savefig("common_mask.pdf")
#f.savefig("common_mask_all_models.pdf")
#f.savefig("common_mask_all_models_w_deserts.pdf")
# -




