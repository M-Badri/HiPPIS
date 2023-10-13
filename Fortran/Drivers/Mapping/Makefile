## Intel compiler configurations
## xMIC-AVX512 is needed for vectorization performance gain on KNL. 
#FC      := ifort   
#F_OPT   :=-O2 -stand f95 -warn none -qopenmp -qopenmp-simd -align array64byte -fp-model strict 
## Use the line below for performance gain on KNL
##F_OPT   :=-O2 -xMIC-AVX512 -warn none -qopenmp -qopenmp-simd -align array64byte -fp-model strict

## gfortran configurations
FC      := gfortran
F_OPT   :=-O2 -Wall -w  -g -pedantic -Wunused-parameter -fcheck=all -std=f95 -ffpe-trap=invalid,zero,overflow -std=f95 -fopenmp -ftree-vectorize 

FFLAGS  := $(F_OPT)
INCLUDE :=-I./

TEST_MAIN:= build
MAIN := main.o nms.o mod_adaptiveInterpolation.o MQSI.o EVAL_BSPLINE.o SPLINE.o REAL_PRECISION.o


all : $(TEST_MAIN)

MQSI.o          : REAL_PRECISION.o lapack.o blas.o
EVAL_BSPLINE.o  : REAL_PRECISION.o
SPLINE.o        : EVAL_BSPLINE.o REAL_PRECISION.o
main.o          : nms.o mod_adaptiveInterpolation.o MQSI.o EVAL_BSPLINE.o SPLINE.o REAL_PRECISION.o
tutorial.o      : mod_adaptiveInterpolation.o
testing.o       : mod_adaptiveInterpolation.o



clean :
	rm -f ../../Src/*.mod ../../Src/*.o
	rm -f ../../Extras/*.mod ../../Extras/*.o
	rm -f *.mod *.out *.o  *.dat main vectorization_results approximations_tables_1d_2d.txt *__* 


$(TEST_MAIN) : $(MAIN) $(TUTORIAL) $(TESTING)
	$(FC) $(FFLAGS) $(MAIN) -o main -lm ../../Extras/lapack.f ../../Extras/blas.f
	mkdir -p mapping_data/data
	ulimit -s unlimited
	./main
	@echo "The stack size has been increased to accommodate 2D examples"
	@echo "The 2D examples take a long time because the solution is evaluated onto a 1000 x 1000 mesh for each example and saved."
	@echo "ran ./main to produce results for the examples in manuscript"
#------------------------------------------------------------------------------
# Default rules
#------------------------------------------------------------------------------
%.o : ../../Extras/%.f90
	$(FC) -c $(FFLAGS) $(INCLUDE) $< -o $@

%.o : ../../Extras/%.F90
	$(FC) -c $(FFLAGS) $(INCLUDE)  $< -o $@

%.o : ../../Extras/%.F
	$(FC) -c $(FFLAGS) $(INCLUDE)  $< -o $@

%.o : ../../Extras/%.f
	$(FC) -c $(FFLAGS) $(INCLUDE)  $< -o $@

%.o : ../../Src/%.F90
	$(FC) -c $(FFLAGS) $(INCLUDE)  $< -o $@

%.o : %.F90
	$(FC) -c $(FFLAGS) $(INCLUDE)  $< -o $@


