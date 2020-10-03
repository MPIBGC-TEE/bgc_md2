{my_python, np}:
    ( 
      with my_python.pkgs; 
          { 
            inherit buildPythonPackage fetchPypi isPyPy numpy cython cftime mpi4py ;
          }
    ) 
    //
    ( 
      with np;
          { 
            inherit openssh libjpeg zlib curl stdenv lib ;
            netcdf=netcdf-mpi; #important for parallel IO version
            hdf5=hdf5-mpi; # important for parallel version
          }
    )
