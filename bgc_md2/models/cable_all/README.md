## Cable versions

This directory contains several versions of cable that can be succesfully built 
with nix.
The one used to produce the data for the transit time and age distributions is: 

## cable_transit_time

If you create a nix-shell environment by typing:

nix-shell

in this directory which will eveluate the expression 

```shell.nix```

This will cause the following executables to be built along with the
dependencies (netcdffortran,mpi...):

```
cable
cable-mpi
``` .

To run cable you will also need appropriate input files.
You can find small script examples in.

```cable_transit_time/example_runs/*.sh```


At the moment all of them use the parallel version cable-mpi
Since this is necessarry for the data we want to use.
