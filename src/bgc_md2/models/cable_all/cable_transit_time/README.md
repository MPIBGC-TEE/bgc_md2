## Build 

```nix-shell ```

Will produce a complete build environment for cable and also the cable-mpi 
executable. By this
```shell.nix``` 
will be evaluated and is a good place to start unravelling 
what happens. It points to the necessary files 
```builder.sh, setup.sh and the sources.

## Run
Some example runs are encoded as scripts in ```example_runs```
Note that the input and output data files are not part of the repository 
but have to be present in the location 
```$cableDataDir```

Which is for now set by the expression in ```shell.nix``` to the location on 
```matagorda``` and ```Antakya``` but can be overriden
by the user to a site specific value.

Note that
the code to fetch the input files for the scripts in ```example_runs/```
is highly dependent on the internal structure
of this ```$cableDataDir``` If you want to reproduce the examples you need 
a copy of this directory with at least all input files present in the same
location. 


## Postprocessing

Notes from China:
For multiyear runs the output of out_cable.nc
has to be renamed (by giving the year in the filename) after every year to prevent the postprocessing script from overiding it
It should be possible to get this out of the namelist file
To run cable for many years you have to restart it again 
For the original look at the file:
My Passport/cable_chris_transit_time/CABLE-run-test-lqy/spinup/run_spinup 
For the cable run as Chris did for the transit time papers we use the  script mknml.bash in cable_chris_transit_time/mm/ on PASSPORT
we have to change the paths in this scripts
instead of /data/lo02b/CABLE-run-tracebility we need  My_Passport/cable_chris_transit_time/CABLE-AUX/offline/gridinfo_NCAR_1.9x2.5_landfrac_revised
veg_params_cable_MK3L_v2_kf.txt
 
 age4_annual.ncl produces fluxes (in this script tha last of before years values are used (364,:,:,:) since they are updated later or eqrlier than other varibles
 see also the pseudo code in  GetFluxFromCABLEoutput.txt in PASSPort/cable_chris_transit_time/CABLE-transit-time-real/archive_TransitTime/scripts (array dimensions do not always match there)
 - in ncl 0 means first index (like python)

