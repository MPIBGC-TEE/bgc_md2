#==============================================================================
TESTFILE=simple_xy_par.nc
# these are all the files we are compiling
#LSRC = $(name).F90
LSRC = simple_xy_par_wr.F90

# this is the executable we are building
FortranProg = $(name) 
PythonProg = "readNcdfPythonParrallel.py"

# from the compiled object files 
OBJS	= ${LSRC:.F90=.o} 

#compiler switches and flags
CINC = -I$(NCMOD)

#suffixes we use
.SUFFIXES:
.SUFFIXES: .F90 .o 
#default rules for these suffixes
.F90.o:
	$(FC) $(CFLAGS) -I$(NCMOD) -c $<

# default target by convention is ``all''
all : $(FortranProg)

#build FortranProg (executable) by linking all objects
$(FortranProg) : $(OBJS)
	$(FC) $(LDFLAGS) -o $@ $(OBJS) $(myLD) 


## test that we can read the file written by the fortran code in python ## the
#oversubscribe flag just allows us to use more processes than we have processors
#(which one would not do on a real cluster for performance reasons) but which is
#fine for testing that we can write and read the same file with a different
#number of processors
testFortran: $(FortranProg) 
	mpirun -n  5 --oversubscribe		$(FortranProg)
	mpirun -n 16 --oversubscribe python3 	$(PythonProg) 

test:
	make clean && make testFortran 

check: test

clean:
	rm -f *.o *.mod $(FortranProg) $(TESTFILE) parallel_test.nc
