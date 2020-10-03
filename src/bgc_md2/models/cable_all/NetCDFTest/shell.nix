let np = import <nixpkgs> {};
in
  np.mkShell {
    buildInputs =
    [
      (
        (import ../my_python_with_packages.nix)
        {
           np = np;
           my_netcdf = np.pkgs.netcdf-mpi;
           my_hdf5 = np.pkgs.hdf5-mpi;
        }
      )
    ]
    ++
    (import ../cable_build_inputs.nix np.pkgs).buildInputs
    ++
    (
      map import [
        ./NCDFParallelIOExample_2D/default.nix
        ./NCDFParallelIOExample_original/default.nix
        ./CableLikeParallelIOExample/default.nix
      ]
    )
    ;
 }
