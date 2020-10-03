{my_python, np}:
  (
    with my_python.pkgs;{ 
      inherit 
        buildPythonPackage 
        fetchPypi 
        numpy 
        pandas 
        pyarrow 
        pytestrunner 
        pytest
        h5py;
      }
  )
  //
  (with np; {inherit lib;})
