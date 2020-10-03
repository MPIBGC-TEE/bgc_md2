let np = import <nixpkgs> {};
stdenv = np.pkgs.stdenv;
in
  stdenv.mkDerivation
  (
    import ../commonShellVars.nix
    //
    {
      name = "NCDFParallelIOExample_2D";
      specificSrcDir = ./src;
    }
  )
