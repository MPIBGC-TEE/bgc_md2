let 
  np = import <nixpkgs> {};
  my_python = np.python37;
in 
  #import ./default.nix ( import ./my_args.nix {my_python=my_python; np=np; })
  import /home/mm/nixpkgs_mm_1/pkgs/development/python-modules/netcdf4/default.nix ( import ./my_args.nix { my_python=my_python; np=np; })
  
