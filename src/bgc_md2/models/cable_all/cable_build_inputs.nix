pkgs:
    let 
    my_netcdf_fortran = pkgs.netcdffortran.override
    (
      with pkgs;{
             netcdf = netcdf-mpi;
             hdf5 = hdf5-mpi;
      }
    );
    # let my_netcdf_fortran = import <nixpkgs/pkgs/development/libraries/netcdf-fortran/default.nix>
    #(
    #  with pkgs;
    #       {
    #           inherit stdenv curl fetchurl ;
    #           # chose the compiler we also use for cable later
    #           gfortran = compiler;

    #           # parallel IO in cable
    #           netcdf = netcdf-mpi;
    #           hdf5 = hdf5-mpi;
    #      }
    #)
    in
    { buildInputs = ( with pkgs;
        [
          gfortran
          findutils
          coreutils
          patchelf
          gnumake
          my_netcdf_fortran
          openmpi
          gdb
          openssh #only needed at runtime
        ]
        );
        my_netcdf_fortran = my_netcdf_fortran;
   }
