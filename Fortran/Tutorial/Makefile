## Intel compiler configurations
## xMIC-AVX512 is need for vectorization performace gain. This means
## The systems supports AVX512 instructions.
#FC      := ifort  
##F_OPT   :=-O3 -xMIC-AVX512 -qopenmp -qopenmp-simd -align array64byte -fp-model strict 
##F_OPT   :=-O2 -warn -stand f90 -qopenmp -align array64byte -fp-model strict 
#F_OPT   :=-O2 -qopenmp -stand f95 -fp-model strict 
#F_PREC  :=-r8

# gfortran configurations
FC      := gfortran
#F_OPT   :=-O3 -Wall -w -g -fcheck=all -fopenmp  -ftree-vectorize -freal-4-real-8 -Wno-argument-mismatch -falign-loops=64 
F_OPT   :=-O2 -Wall -w -g -fcheck=all -std=f95 -ffpe-trap=invalid,zero,overflow -std=f95 -fopenmp 

FFLAGS  := $(F_OPT)
INCLUDE :=-I./

TEST_MAIN:= build
TUTORIAL := tutorial.o mod_adaptiveInterpolation.o
TESTING := testing.o mod_adaptiveInterpolation.o


all : $(TEST_MAIN)

tutorial.o      : mod_adaptiveInterpolation.o
testing.o       : mod_adaptiveInterpolation.o



clean :
	rm -f $(OBJS) *.mod *.out *.o  *.dat testing tutorial 


$(TEST_MAIN) : $(MAIN) $(TUTORIAL) $(TESTING)
	$(FC) $(FFLAGS) $(TUTORIAL) -o tutorial
	$(FC) $(FFLAGS) $(TESTING) -o testing
	ulimit -s unlimited
	@echo 
	@echo "run ./tutorial to execute tutorial examples"
	@echo 
	@echo "run ./testing to execute the test examples"
	@echo 
#------------------------------------------------------------------------------
# Default rules
#------------------------------------------------------------------------------
%.o : %.F90
	$(FC) -c $(FFLAGS) $(INCLUDE)  $< -o $@

