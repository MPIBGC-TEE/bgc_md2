let np = import <nixpkgs> {};
stdenv = np.pkgs.stdenv;
in
  stdenv.mkDerivation
  (
    import ../commonShellVars.nix
    //
    {
      name = "CableLikeParallelIOExample";
      specificSrcDir=./src;
    }
  )
