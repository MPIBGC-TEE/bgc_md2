#!/usr/bin/env python
import numpy as np
import general_helpers as gh
import matplotlib.pyplot as plt
from pathlib import Path


model_folders = [
    'yz_jules',
    'cj_isam',
    'Aneesh_SDGVM',
    'kv_ft_dlem',
    'kv_visit2',
    'jon_yib',
    'bian_ibis2'
]


def mask(mf):
    print("###########################################")
    print("#######   {}   ######".format(mf))
    msh = gh.msh(mf)
    return msh.spatial_mask(Path(gh.confDict(mf)['dataPath']))

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
common_mask = gh.combine_masks(masks)
# common_mask=gh.project_2(source=masks[0],target=masks[1])
# common_mask=gh.project_2(source=masks[0],target=ut)
common_mask.write_netCDF4(Path("common_mask.nc"))
ax = f.add_subplot(n+1, 1, n+1)
common_mask.plot_dots(ax)


f.savefig("common_mask.pdf")
