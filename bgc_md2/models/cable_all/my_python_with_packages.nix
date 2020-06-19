{np,my_netcdf,my_hdf5}:
  let
    my_python = np.python38;
    my_args =
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
        }
    )
    //
    {
      netcdf=my_netcdf;
      hdf5=my_hdf5;
    }
    ;
    my_netcdf4_python = import ./my_netcdf4_python/default.nix my_args ;
  in
    (my_python.withPackages ((ps: [
      ps.mpi4py
      ps.numpy
      ps.sympy
      ps.bootstrapped-pip
      ps.jupyterlab
      ps.ipyparallel # jupyter serverextension enable --py ipyparallel
      #ps.netcdf4
    ]
    ++ [my_netcdf4_python]
    )))

