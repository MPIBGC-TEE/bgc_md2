let
    np = import <nixpkgs> {};
    pkgs=np.pkgs;
    my_netcdf=pkgs.netcdf-mpi;
    my_hdf5=pkgs.hdf5-mpi;
    cable_build_inputs = (import ../cable_build_inputs.nix  pkgs);

    # add some python tools
    my_python = np.python38;

    my_python_with_packages = import ../my_python_with_packages.nix
    {
      np = np;
      my_netcdf = my_netcdf;
      my_hdf5 = my_hdf5;
    };


in {
    buildInputs =
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
    FC="mpif90";
    NCMOD="${cable_build_inputs.my_netcdf_fortran}/include";
    CFLAGS="-x f95-cpp-input -g ";
    #CFLAGS="-x f95-cpp-input";
    myLD="-lnetcdff";
    commonSrcDir =./commonSrc;
    #name="simple_xy_par_wr";
    unpackPhase = ''
      echo "################ unpack phase"
      cp $specificSrcDir/* .
      cp $commonSrcDir/Makefile .
      cp $commonSrcDir/readNcdfPythonParrallel.py .
      ls
    '';

    buildPhase=''
      echo "################ build phase"
      make -f Makefile
    '';

    doCheck = true;
    checkPhase=''
      echo '################ check phase'; make -f Makefile test
    '';

    installPhase=''
    mkdir -p $out/bin/
    cp $name $out/bin/
    '';
}
