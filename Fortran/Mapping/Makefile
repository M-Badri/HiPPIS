## Intel compiler configurations
## xMIC-AVX512 is need for vectorization performace gain. This means
## The systems supports AVX512 instructions.
#FC      := ifort  
#F_OPT   :=-O3 -xMIC-AVX512 -qopenmp -qopenmp-simd -align array64byte -fp-model strict 
#F_PREC  :=-r8

# gfortran configurations
FC      := gfortran
F_OPT   :=-O3 -Wall -w -g -fcheck=all -fopenmp  -ftree-vectorize -freal-4-real-8 -Wno-argument-mismatch -falign-loops=64 
F_PREC  :=-fdefault-real-8

FFLAGS  := $(F_OPT)
INCLUDE :=-I./

TEST_MAIN:= build
MAIN := main.o nms.o mod_adaptiveInterpolation.o
TUTORIAL := tutorial.o mod_adaptiveInterpolation.o
TESTING := testing.o mod_adaptiveInterpolation.o


all : $(TEST_MAIN)


main.o        : nms.o mod_adaptiveInterpolation.o
tutorial.o    : mod_adaptiveInterpolation.o
testing.o     : mod_adaptiveInterpolation.o


clean :
	rm -f $(OBJS) *.mod *.out *.o  *.dat


$(TEST_MAIN) : $(MAIN) $(TUTORIAL) $(TESTING)
	$(FC) $(FFLAGS) $(MAIN) -o main
	$(FC) $(FFLAGS) $(TUTORIAL) -o tutorial
	$(FC) $(FFLAGS) $(TESTING) -o testing
	mkdir -p mapping_data/data
	ulimit -s unlimited
	@echo "the stack size has been increased to accomodate 2D examples"
	@echo "The 2D examples take long time because the solution is evaluated onto a 1000 x 1000 mesh for each example and saved."
	./main
	@echo "ran ./main to produce results for the examples in manuscrip"
	@echo 
	@echo "run ./tutorial to execute tutorial examples"
	@echo 
	@echo "run ./testing to execute the test examples"
	@echo 
#------------------------------------------------------------------------------
# Default rules
#------------------------------------------------------------------------------
%.o : %.f90
	$(FC) -c $(FFLAGS) $(INCLUDE) $< -o $@

%.o : %.F90
	$(FC) -c $(FFLAGS) $(INCLUDE)  $< -o $@
