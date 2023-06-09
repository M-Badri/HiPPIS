## Intel compiler configurations
## xMIC-AVX512 is need for vectorization performace gain. This means
## The systems supports AVX512 instructions.
#FC      := ifort  
#F_OPT   :=-O2 -warn -stand f95 -qopenmp -qopenmp-simd -align array64byte -fp-model strict 
#F_PREC  :=-r8

# gfortran configurations
FC      := gfortran
F_OPT   :=-O2 -Wall -w -g -fcheck=all -std=f95 -ffpe-trap=invalid,zero,overflow -fopenmp 

FFLAGS  := $(F_OPT)
INCLUDE :=-I./

TEST_MAIN:= build
TUTORIAL := main.o mod_adaptiveInterpolation.o
TESTING := testing.o mod_adaptiveInterpolation.o


all : $(TEST_MAIN)

main.o      : mod_adaptiveInterpolation.o
testing.o       : mod_adaptiveInterpolation.o



clean :
	rm -f $(OBJS) *.mod *.out *.o  *.dat testing main 


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

