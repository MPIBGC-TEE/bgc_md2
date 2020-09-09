let 
  np = import <nixpkgs> {};
  my_python = np.python38;
  my_awkward= import /home/mm/nixpkgs_mm/pkgs/development/python-modules/awkward ( import ./test_args.nix {my_python=np.python37; np=np; });
  
in np.stdenv.mkDerivation {
    name ="test";
    buildInputs = with np.pkgs; [ 
      (my_python.withPackages ((ps: 
      [
        ps.ipython
      #  ps.bootstrapped-pip 
      #  ps.jupyterlab
      #  ps.ipyparallel # necessarry
      #  ps.dask
      #  ps.xarray
      ##  ps.netcdf4 alternative to 
      ] 
      ++ 
      [my_awkward]
    )))
    ];
}
