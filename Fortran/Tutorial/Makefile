## Intel compiler configurations
## on -xMIC-AVX512 flag is needed for vectorization performace gain.
#FC      := ifort  
#F_OPT   :=-O2 -warn -stand f95 -qopenmp -qopenmp-simd -align array64byte -fp-model strict 


## gfortran configurations
#FC      := gfortran
#F_OPT   :=-O2 -Wall -w -g -fcheck=all -std=f95 -ffpe-trap=invalid,zero,overflow -fopenmp 

# NAG configuration
FC      := /home/tajo10/usr/local/NAG_Fortran/bin/nagfor
F_OPT   := -gline #-g -u  -nan -C=undefined -C=dangling -kind=byte -dusty
#F_OPT   :=-O3 -g -g90 -openmp=4.5 -u -nan -C=all -C=undefined -C=dangling -kind=byte -dusty

FFLAGS  := $(F_OPT)
INCLUDE :=-I./

TEST_MAIN:= build
TUTORIAL := main.o mod_adaptiveInterpolation.o
TESTING := testing.o mod_adaptiveInterpolation.o


all : $(TEST_MAIN)

main.o      : mod_adaptiveInterpolation.o
testing.o   : mod_adaptiveInterpolation.o



clean :
	rm -f $(OBJS) *.mod *.out *.o  *.dat *.mex* testing main 


$(TEST_MAIN) : $(MAIN) $(TUTORIAL) $(TESTING)
	$(FC) $(FFLAGS) $(TUTORIAL) -o main
	$(FC) $(FFLAGS) $(TESTING) -o testing
	@echo 
	@echo "run ./main to execute tutorial examples in main.F90"
	@echo 
	@echo "run ./testing to execute the test examples"
	@echo 
#------------------------------------------------------------------------------
# Default rules
#------------------------------------------------------------------------------
%.o : %.F90
	$(FC) -c $(FFLAGS) $(INCLUDE)  $< -o $@

