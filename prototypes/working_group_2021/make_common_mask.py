from pathlib import Path
import numpy as np
import general_helpers as gh
import matplotlib.pyplot as plt


model_folders=['kv_visit2','cj_isam']#'Aneesh_SDGVM','jon_yib' ,'cable-pop','kv_visit2','cj_isam','yz_jules','kv_ft_dlem']
def mask(mf):
    msh = gh.msh(mf)
    return msh.spatial_mask(Path(gh.confDict(mf)['dataPath']))
    
masks=list(
        map(
            mask,
            model_folders
        )
)
ma=np.zeros(
    shape=(360,720),
    #dtype=np.bool_
)
ut = gh.CoordMask(
        index_mask=ma,
        tr=gh.globalMaskTransformers(ma)
)
#masks.append(ut)
#############################################################################################
#plotting
f=plt.figure(figsize=(10,20))
n=len(masks)
for i in range(n):
    m=masks[i]
    ax=f.add_subplot(n+1,1,i+1)
    m.plot_dots(ax)

#At the end of the list add the  ultimate target mask (empty)
common_mask=gh.combine_masks(masks)
#common_mask=gh.project_2(source=masks[0],target=masks[1])
#common_mask=gh.project_2(source=masks[0],target=ut)
#common_mask.write_netCDF4(Path("common_mask.nc"))
ax=f.add_subplot(n+1,1,n+1)
common_mask.plot_dots(ax)


f.savefig("common_mask.pdf")
