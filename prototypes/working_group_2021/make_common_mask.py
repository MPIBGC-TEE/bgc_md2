from pathlib import Path
import numpy as np
import general_helpers as gh


model_folders=['kv_visit2','cj_isam']#,'jon_yib' ,'cable-pop','kv_visit2','cj_isam','yz_jules','kv_ft_dlem']
def mask(mf):
    msh = gh.msh(mf)
    return msh.spatial_mask(Path(gh.confDict(mf)['dataPath']))
    
m1=mask('kv_visit2')
m2=mask('cj_isam')
print(m1.index_mask.mean())
print(m2.index_mask.mean())
masks=[m1,m2]
#masks=list(
#        map(
#            mask,
#            model_folders
#        )
#)
ma=np.zeros(
    shape=(360,720),
    #dtype=np.bool_
)
ut = gh.CoordMask(
        index_mask=ma,
        tr=gh.globalMaskTransformers(ma)
)
#masks.append(ut)
for m in masks:
    print(m.index_mask)
    print(m.index_mask.mean())
#At the end of the list add the  ultimate target mask (empty)
#common_mask=gh.combine_masks(masks+[ut])
common_mask=gh.project_2(source=masks[0],target=masks[1])
#common_mask=gh.project_2(source=masks[0],target=ut)
common_mask.write_netCDF4(Path("common_mask.nc"))
