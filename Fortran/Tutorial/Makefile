# Intel compiler configurations
# on -xMIC-AVX512 flag is needed for vectorization performace gain on KNL.
FC      := ifort  
#F_OPT   :=-O2 -warn -qopenmp -qopenmp-simd -align array64byte -fp-model strict 
F_OPT   :=-O2 -warn -stand f90 -qopenmp -qopenmp-simd -align array64byte -fp-model strict 


## gfortran configurations
#FC      := gfortran
#F_OPT   :=-O2 -Wall -w -g -fcheck=all -std=f95 -ffpe-trap=invalid,zero,overflow -fopenmp 

FFLAGS  := $(F_OPT)
INCLUDE :=-I./

TEST_MAIN:= build
TUTORIAL := main.o mod_adaptiveInterpolation.o
TESTING := testing.o mod_adaptiveInterpolation.o


all : $(TEST_MAIN)

main.o      : mod_adaptiveInterpolation.o
testing.o   : mod_adaptiveInterpolation.o



clean :
	rm -f $(OBJS) *.mod *.out *.o  *.dat *.mex* testing main *__* 


$(TEST_MAIN) : $(MAIN) $(TUTORIAL) $(TESTING)
	$(FC) $(FFLAGS) $(TUTORIAL) -o main
	$(FC) $(FFLAGS) $(TESTING) -o testing
	./main
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

