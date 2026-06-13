let 
  np = import <nixpkgs> {};
in 
  import /home/mm/nixpkgs_mm/pkgs/development/python-modules/awkward ( import ./test_args.nix {my_python=np.python37; np=np; })
