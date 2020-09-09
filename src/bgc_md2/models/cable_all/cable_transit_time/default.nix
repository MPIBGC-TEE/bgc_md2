let
    np = import <nixpkgs> {};
    pkgs=np.pkgs;
    cable_build_inputs = (import ../cable_build_inputs.nix  pkgs);

in builtins.derivation {
    builder = "${pkgs.bash}/bin/bash";
    args= [ ./builder.sh ];
    system =builtins.currentSystem;
    name = "cable";
    #src = null;
    src = np.lib.cleanSource ./CABLE-SRC;
    setup =./setup.sh;
    buildInputs = cable_build_inputs.buildInputs;
    ncd= cable_build_inputs.my_netcdf_fortran;
}

