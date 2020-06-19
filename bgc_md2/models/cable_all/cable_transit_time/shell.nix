let
    np = import <nixpkgs> {};
    pkgs=np.pkgs;
    my_netcdf=pkgs.netcdf-mpi;
    my_hdf5=pkgs.hdf5-mpi;
    cable_build_inputs = (import ../cable_build_inputs.nix  pkgs);

    cable = import ./default.nix;

    # add some python tools
    my_python = np.python38;

    my_python_with_packages = import ../my_python_with_packages.nix 
    {
      np = np; 
      my_netcdf = my_netcdf;
      my_hdf5 = my_hdf5;
    };


in np.mkShell {
    buildInputs =
      [cable]
      ++
      cable_build_inputs.buildInputs
      # add some extra tooling
      ++
      [my_python_with_packages my_netcdf]
      ++
      ( with np.pkgs;
        [
          openssh which gnutar gzip bzip2 gcc binutils-unwrapped
          gawk gnused gnugrep
          nodejs #necessary for jupyter lab vi extension
        ]
      );
    #set some variables for convinience 
    cableDataDir="/home/data/cable-data/example_runs"; #site specific
    np=32 ; #seems fastest
    cable_exe="cable-mpi";
}
