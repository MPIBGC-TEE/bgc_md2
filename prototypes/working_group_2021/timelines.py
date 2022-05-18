# ---
# jupyter:
#   jupytext:
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

import netCDF4 as nc
import datetime as dt
import general_helpers as gh
from pathlib import Path
#############################################################################
# SDGVM
#############################################################################
time=nc.Dataset(
    Path(gh.confDict("Aneesh_SDGVM")["dataPath"]).joinpath("SDGVM_S2_npp.nc")
).variables["time"]
print(time)
start_date=dt.date(1901,1,1)
# interpretation 1, believe:
# -the unit as hours since 1901,1,1
delta_1=dt.timedelta(hours=float(time[-1].data))
end_date_1=start_date+delta_1
print("\n\nEnd date according to unit hours since 1900-01-01: {}".format(start_date+delta_1))
# interpretation 2, believe
# -the dimension of the timevariable:number of months=1440, 
# -the number of days per month=30
delta_2=dt.timedelta(days=time.shape[0]*30)
print("\n\nEnd date according to number of 30 day months since 1900-01-01: {}".format(start_date+delta_2))






