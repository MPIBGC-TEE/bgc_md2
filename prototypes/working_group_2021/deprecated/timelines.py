# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.13.8
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
reference_date=dt.date(1900,1,1)
first_value_in_days=time[0]/24
print("\nFirst value in days: " + str(first_value_in_days))
start_date=reference_date+dt.timedelta(first_value_in_days)
print("Start date: "+str(start_date))
last_value_in_days=time[-1]/24
print("\nLast value in days: " + str(last_value_in_days))
end_date_real_months=reference_date+dt.timedelta(last_value_in_days)
print("End date (with real months): "+str(end_date_real_months))
last_date_in_days=reference_date.year*360+reference_date.month*30+reference_date.day+last_value_in_days
end_year_30d=last_date_in_days//360
end_month_30d=(last_date_in_days%360)//30
end_day_30d=(last_date_in_days%360)%30
print("End date (with 30day months): "+str(end_year_30d)+"-"+str(end_month_30d)+"-"+str(end_day_30d))
print ("2020-00-16 actually means 2019-12-16 - this is exactly what Panoply shows as the date of last measurement")

# interpretation 1, believe:
# -the unit as hours since 1901,1,16
delta_1=dt.timedelta(hours=float(time[-1].data)+15*24) # 15 days because of the last recorded time being mid month
end_date_1=start_date+delta_1
print("\n\nEnd date according to unit hours since 1900-01-16, disregarding 12 months/year: {}".format(start_date+delta_1))

# interpretation 2, believe:
# -the dimension of the timevariable:number of months=1440, 
# -the number of days per month=30
delta_2=dt.timedelta(days=time.shape[0]*30)
print("\n\nEnd date according to number of 30 day months since 1900-01-16, disregarding 12 months/year: {}".format(start_date+delta_2))
# note that this is actually compatible with the former result 

# interpretation 3, believe
# -the dimension of the timevariable:number of months=1440, 
# -the number of months/year=12
print("\n\nEnd date according to 1440 months and 12 month years since 1900-01-16: {}".format(start_date.year+1440/12))

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
# -

# now import the NPP data and explore how it actually behaves
import sys
sys.path.insert(0,'Aneesh_sdgvm') # necessary to import general_helpers
import model_specific_helpers_2 as msh
from pathlib import Path
import json
with Path('Aneesh_sdgvm/config.json').open(mode='r') as f:
    conf_dict=json.load(f) 
svs,dvs=msh.get_global_mean_vars(dataPath=Path(conf_dict["dataPath"]))

# now we plot distribution of NPP values over 12 months
import matplotlib.pyplot as plt
import numpy as np
# 1st plot is the beginning of the simulation
plt.plot(
    dvs.npp[0:12])
plt.xticks(np.arange(0, 12, 1))
plt.grid()

# 2nd plot is the middle of the simulation
plt.plot(
    dvs.npp[len(dvs.npp)//2-12:len(dvs.npp)//2])
plt.xticks(np.arange(0, 12, 1))
plt.grid()

# 3rd plot is the end of the simulation
plt.plot(
    dvs.npp[len(dvs.npp)-12:len(dvs.npp)])
plt.xticks(np.arange(0, 12, 1))
plt.grid()

# +
# we see that the npp follows the months throughout the simulation, so interpretation 3 is correct
# -


