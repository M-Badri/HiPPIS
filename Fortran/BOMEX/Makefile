## intel compiler configurations
FC      := ifort
F_OPT   :=-O2 -fp-model  strict -stand f95 -ftz -fpe0 -traceback -qopenmp-simd -align array64byte  

## gfortran compiler configurations
#FC     := gfortran
#F_OPT   :=-O2 -g -Wall -fcheck=all -std=f95 -fopenmp -fbacktrace -w -ffpe-trap=invalid,zero,overflow  


FFLAGS  := $(F_OPT)
INCLUDE :=-I./

MAIN    := main

OBJS_BASIC := machine.o physcons.o #mod_advection.o
MQS_OBJS := MQSI.o EVAL_BSPLINE.o SPLINE.o REAL_PRECISION.o
OBJS := main.o mod_bomex.o mod_adaptiveInterpolation.o nms.o mfpbl.o moninedmf.o scm_physics.o kessler.o $(MQS_OBJS) $(OBJS_BASIC)
all : $(MAIN)

mod_bomex.o    : mod_adaptiveInterpolation.o $(OBJS_BASIC)
moninedmf.o    : mfpbl.o $(OBJS_BASIC)
mfpbl.o        : $(OBJS_BASIC)
kessler.o      : $(OBJS_BASIC)
scm_physics.o  : mod_adaptiveInterpolation.o kessler.o moninedmf.o $(OBJS_BASIC)
MQSI.o          : REAL_PRECISION.o lapack.o blas.o
EVAL_BSPLINE.o  : REAL_PRECISION.o
SPLINE.o        : EVAL_BSPLINE.o REAL_PRECISION.o
main.o         : mod_bomex.o moninedmf.o scm_physics.o mod_adaptiveInterpolation.o $(MQS_OBJS) $(OBJS_BASIC)


clean :
	rm -f $(OBJS) *.mod main *__* *.o

$(MAIN) : $(OBJS)
	$(FC) $(LDFLAGS) $(OBJS) -o $(MAIN) -lm lapack.f blas.f
	mkdir -p bomex_data
	echo 'The Bomex example takes a long time. Run ./main'
	./$(MAIN)
        

#------------------------------------------------------------------------------
# Default rules
#------------------------------------------------------------------------------
%.o : %.f90
	$(FC) -c $(FFLAGS) $(INCLUDE) $< -o $@

%.o : %.F90
	$(FC) -c $(FFLAGS) $(INCLUDE) $(CPPFLAGS) $< -o $@
%.o : %.f
	$(FC) -c $(FFLAGS) $(INCLUDE) $(CPPFLAGS) $< -o $@
