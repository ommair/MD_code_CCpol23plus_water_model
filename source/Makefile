# This is Makefile for creating executables

# This variable is for fortran compiler
FC = mpif90 #mpiifort #mpif90 #gfortran

# This variable will give compiler flags
FFLAGS = -O3 #-Ofast 
FCFLAGS = -c 
FLAGS = #-fbounds-check

OBJS = common_arrays.o \
       utility.o       \
       molecular_sites.o \
       ewald.o          \
       vdw_lj.o            \
       vdw_sapt5s.o     \
       vdw_ccpol5s.o    \
       vdw_ccpol8s.o    \
       vdw_ccpol8s_omo.o    \
       vdw_ccpol23plus.o     \
       nb_induction_model.o  \
       threebody_potential.o  \
       statistics.o     \
       initialization.o    \
       evolve.o        \

EXECUTE = execute

$(EXECUTE): md.f90 $(OBJS)
	$(FC) $(FFLAGS) -o $(EXECUTE) md.f90 $(OBJS)
 
# These lines produces the .mod and .o files
common_arrays.o: common_arrays.f90
	$(FC) $(FLAGS) $(FFLAGS) $(FCFLAGS) common_arrays.f90
utility.o: utility.f90
	$(FC) $(FLAGS) $(FFLAGS) $(FCFLAGS) utility.f90
molecular_sites.o: molecular_sites.f90
	$(FC) $(FLAGS) $(FFLAGS) $(FCFLAGS) molecular_sites.f90
initialization.o: initialization.f90
	$(FC) $(FLAGS) $(FFLAGS) $(FCFLAGS) initialization.f90
vdw_lj.o: vdw_lj.f90
	$(FC) $(FLAGS) $(FFLAGS) $(FCFLAGS) vdw_lj.f90
vdw_sapt5s.o: vdw_sapt5s.f90 
	$(FC) $(FLAGS) $(FFLAGS) $(FCFLAGS) vdw_sapt5s.f90
vdw_ccpol5s.o: vdw_ccpol5s.f90
	$(FC) $(FLAGS) $(FFLAGS) $(FCFLAGS) vdw_ccpol5s.f90
vdw_ccpol8s.o: vdw_ccpol8s.f90
	$(FC) $(FLAGS) $(FFLAGS) $(FCFLAGS) vdw_ccpol8s.f90
vdw_ccpol8s_omo.o: vdw_ccpol8s_omo.f90
	$(FC) $(FLAGS) $(FFLAGS) $(FCFLAGS) vdw_ccpol8s_omo.f90
vdw_ccpol23plus.o: vdw_ccpol23plus.f90
	$(FC) $(FLAGS) $(FFLAGS) $(FCFLAGS) vdw_ccpol23plus.f90 
evolve.o: evolve.f90 
	$(FC) $(FLAGS) $(FFLAGS) $(FCFLAGS) evolve.f90
statistics.o: statistics.f90
	$(FC) $(FLAGS) $(FFLAGS) $(FCFLAGS) statistics.f90
ewald.o: ewald.f90
	$(FC) $(FLAGS) $(FFLAGS) $(FCFLAGS) ewald.f90
nb_induction_model.o: nb_induction_model.f90
	$(FC) $(FLAGS) $(FFLAGS) $(FCFLAGS) nb_induction_model.f90
threebody_potential.o: threebody_potential.f90
	$(FC) $(FLAGS) $(FFLAGS) $(FCFLAGS) threebody_potential.f90
#core_shell.o:core_shell.f90
#	$(FC) $(FLAGS) $(FFLAGS) $(FCFLAGS) core_shell.f90

# This line is for clean up
clean: 
	rm -rf *.o *.mod $(OBJS) *~ $(EXECUTE) OUTPUT *.out *.xyz *.err fort.401 numerical_vs_analytical_forces numerical_vs_analytical_pressure	
