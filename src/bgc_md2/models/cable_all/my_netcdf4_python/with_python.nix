let 
  np = import <nixpkgs> {};
  my_python = np.python38;
  #my_netcdf4_python=import ./default.nix ( import ./my_args.nix {
  my_netcdf4_python=import /home/mm/nixpkgs_mm/pkgs/development/python-modules/netcdf4/default.nix ( import ./my_args.nix {      
    my_python=my_python; 
    np=np; 
  });
  
in np.stdenv.mkDerivation {
    name ="test";
    buildInputs = with np.pkgs; [ 
      (my_python.withPackages ((ps: 
      [
        ps.mpi4py 
        ps.numpy 
        ps.sympy 
        ps.bootstrapped-pip 
        ps.jupyterlab
        ps.ipyparallel # necessarry
        ps.dask
        ps.xarray
      #  ps.netcdf4 alternative to 
      ] 
      ++ 
      [my_netcdf4_python]
    )))
    ];
}
