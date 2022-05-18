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

# +
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
delta_1=dt.timedelta(hours=float(time[-1].data)+15*24) # 15 days because of the last recorded time being mid month
end_date_1=start_date+delta_1
print("\n\nEnd date according to unit hours since 1900-01-01, disregarding 12 months/year: {}".format(start_date+delta_1))

# interpretation 2, believe:
# -the dimension of the timevariable:number of months=1440, 
# -the number of days per month=30
delta_2=dt.timedelta(days=time.shape[0]*30)
print("\n\nEnd date according to number of 30 day months since 1900-01-01, disregarding 12 months/year: {}".format(start_date+delta_2))
# note that this is actually compatible with the former result 

# interpretation 3, believe
# -the dimension of the timevariable:number of months=1440, 
# -the number of months/year=12
print("\n\nEnd date according to 1440 months and 12 month years since 1900-01-01: {}".format(start_date.year+1440/12))

# The consistency of interpretations 1 and 2 suggests that SVGDM runs on physical units and 
# considers a "Month" to be no more than a 30 day time step. 
# - For the Drivers this means that we have a model specific 
#   day2month function: lambda day:int(day/30) which interprets 'monthly' data as the data belonging to a particular
#   30 day period....
#   A consequence is that the seasons would shift which should be visible in the npp data near the end of the simulation
#   E.g. The last 30 day month of SVGDM is between is +-15days of 03-28-2019 on the modern calender 
#   which should make the data look different from the January data. 
# - For the plotting along a calendarial timeline we would have to count the time from the beginning in a physical unit
#   (iterator already implements unit days) and then compute the date (with pythons datetime module).
