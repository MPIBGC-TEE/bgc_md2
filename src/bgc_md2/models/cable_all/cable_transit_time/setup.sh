# This script is called by nix
# and sets up the evironment
unset PATH
for p in $baseInputs $buildInputs; do
  if [ -d $p/bin ]; then
    export PATH="$p/bin${PATH:+:}$PATH"
  fi
	#  if [ -d $p/include ]; then
	#    export NIX_CFLAGS_COMPILE="-I $p/include${NIX_CFLAGS_COMPILE:+ }$NIX_CFLAGS_COMPILE"
	#  fi
	#  if [ -d $p/lib ]; then
	#    export NIX_LDFLAGS="-rpath $p/lib -L $p/lib${NIX_LDFLAGS:+ }$NIX_LDFLAGS"
	#  fi
done
# we export some environmentvariables needed by our version "build_nix.ksh"# of build.ksh  
# we could also change the makefile to directly incorporate them
export NCDIR="$ncd/lib${NCDIR:+ }$NCDIR"
export NCMOD="$ncd/include${NCMOD:+ }$NCMOD"

# define some functions 
function fetchPhase(){
  echo "########################\n \
  copying the following files:"
  echo $src
}
function prepareBuild(){
  echo "########################\n \
  prepareBuild:"
  declare -xp
  #export CFLAGS='-g -O0 -x f95-cpp-input -Wall'
  #export CFLAGS='-g -O0 -x f95-cpp-input' debug
  export CFLAGS='-x f95-cpp-input'
  export LD='-lnetcdff'
  export LDFLAGS="-L ${NCDIR} -O2"
  echo "############# tmpDir=${tmpDir}"
  mkdir -p $tmpDir
  cp -p "${src}/core/biogeophys"/*.F90 $tmpDir
  cp -p "${src}/offline"/*.F90 $tmpDir
  cp -p "${src}/core/biogeochem"/*.F90 $tmpDir
}
function buildSerial(){
  echo "########################\n \
  build serial executable \n \
  #################################\n\n"
  tmpDir="${out}/tmpSerial"
  prepareBuild
  export FC=gfortran
  cp -p ${src}/offline/Makefile_offline  ${tmpDir}
  cd ${tmpDir}
  make -f Makefile_offline
}
function buildParallel(){
  echo "########################\n \
  build parallel executable \n \
  #################################\n\n"
  tmpDir="${out}/tmpParallel"
  prepareBuild
  export FC=mpif90
  cp -p ${src}/offline/Makefile_mpi ${tmpDir}
  cd ${tmpDir}
  make -f Makefile_mpi
}

function installSerial(){
	echo "########################\n \
    install serial executable \n \
    #################################\n\n"
  echo $out
  td=${out}/bin
  mkdir -p ${td}
	mv cable "${td}/"
}
function installParallel(){
	echo "########################\n \
    install parallel executable \n \
    #################################\n\n"
  echo $out
  td=${out}/bin
  mkdir -p ${td}
	mv cable-mpi "${td}/"
}
function fixupPhase(){
	#clean_build
	find $out -type f -exec patchelf --shrink-rpath '{}' \;-exec strip '{}' \; 2>/dev/null
}
function genericBuild(){
	fetchPhase
	buildParallel
	installParallel
	fixupPhase
	buildSerial
	installSerial
	fixupPhase
} 
