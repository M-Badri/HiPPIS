FC     := gfortran 
CC     := gcc
LD     := gfortran
#F_OPT   :=-O2 -g -qopenmp -qopenmp-simd -xCORE-AVX2 -xMIC-AVX512 -traceback -fp-model strict -ftz -fpe0 -align array64byte
#F_OPT   :=-O2 -qopenmp -qno-openmp-simd -no-vec -check bounds -xMIC-AVX512 -traceback -fp-model strict -ftz -fpe0 -align array64byte
#F_OPT   :=-O2 -qopenmp -check bounds -traceback -fp-model strict -ftz -fpe0 -align array64byte
F_OPT   :=-O2 -fopenmp -fbounds-check -g -fcheck=all#-fp-model strict -ftz -fpe0 -align array64byte
F_PREC  :=-r8

FFLAGS  := $(F_OPT)
CPPFLAGS := -fpp
INCLUDE :=-I./

TEST_INTERP := test2
#TEST_SIAM := test
OBJS_BASIC := machine.o 
INTERP := testInterpolation.o nms.o mod_legendre.o mod_adaptiveInterpolation.o $(OBJS_BASIC)
INTERP2 := testInterpolation.f90 nms.f90 mod_legendre.f90 mod_adaptiveInterpolation.f90 machine.f90
SIAM := mainSIAMPaper.o nms.o mod_legendre.o mod_adaptiveInterpolation.o $(OBJS_BASIC)

all : $(TEST_SIAM) $(TEST_INTERP)


testInterpolation.o        : nms.o machine.o mod_legendre.o mod_adaptiveInterpolation.o
mainSIAMPaper.o        : nms.o machine.o mod_legendre.o mod_adaptiveInterpolation.o
mod_legendre.o : $(OBJS_BASIC)


clean :
	rm -f $(OBJS) *.mod a.out siam.out *.o testFunc testInterp a.out *.dat


#$(TEST_SIAM) : $(SIAM)
#	$(LD) $(LDFLAGS) $(SIAM) -o siam.out

$(TEST_INTERP) : $(INTERP)
	$(LD) $(FFLAGS) $(INTERP) -o a.out
	#./a.out &> out
	#mv T* mass_conservation/data/
	#mv u* mass_conservation/data/
	#mv v* mass_conservation/data/
	#mv qv* mass_conservation/data/
	#mv qc* mass_conservation/data/
	#mv pres* mass_conservation/data/
	#mv Heaviside* math_paper_data
	#mv peak* math_paper_data
	#mv Runge* math_paper_data
	#mv eps* math_paper_data
	#mv GelbT* math_paper_data
	#mv Tanh* math_paper_data
	#mv Square* math_paper_data
	#mv T1* math_paper_data
	#mv T2* math_paper_data

#$(TEST_INTERP) : 
#	$(LD) $(FFLAGS) $(INTERP2) -o a.out
#------------------------------------------------------------------------------
# Default rules
#------------------------------------------------------------------------------
%.o : %.f90
	$(FC) -c $(FFLAGS) $(INCLUDE) $< -o $@

%.o : %.F90
	$(FC) -c $(FFLAGS) $(INCLUDE)  $< -o $@
	#$(FC) -c $(CPPFLAGS) $(FFLAGS) $(INCLUDE)  $< -o $@
